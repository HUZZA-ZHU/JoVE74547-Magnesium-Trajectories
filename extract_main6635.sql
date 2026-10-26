-- Reference reconstruction of the FINAL approved cohort algorithm.
-- Engine: DuckDB 1.5.6. Inputs are read-only MIMIC-IV 3.1 raw tables and
-- the official mimic-code Sepsis-3 derived table. All created objects are TEMP.
-- This does not reuse the earlier >=3-raw-record / 6,637-person definition.
-- No trajectory classification or model fitting occurs here.

CREATE OR REPLACE TEMP TABLE final_first_icu AS
SELECT i.subject_id, i.hadm_id, i.stay_id, i.intime, i.outtime,
       p.gender, p.anchor_age + year(i.intime) - p.anchor_year AS age,
       p.dod, a.deathtime, a.dischtime, a.hospital_expire_flag
FROM mimiciv_icu.icustays AS i
JOIN mimiciv_hosp.patients AS p USING (subject_id)
JOIN mimiciv_hosp.admissions AS a USING (subject_id, hadm_id)
QUALIFY row_number() OVER (
  PARTITION BY i.subject_id ORDER BY i.intime, i.stay_id
) = 1;

CREATE OR REPLACE TEMP TABLE final_sepsis_mother AS
SELECT f.*, s.antibiotic_time, s.culture_time,
       s.suspected_infection_time, s.sofa_time,
       s.sofa_score AS sepsis_sofa_score, s.sepsis3,
       f.intime + INTERVAL '120' HOUR AS landmark_time,
       f.intime + INTERVAL '28' DAY AS endpoint_time
FROM final_first_icu AS f
JOIN mimiciv_derived.sepsis3 AS s USING (subject_id, stay_id)
WHERE s.sepsis3 IS TRUE AND f.age BETWEEN 18 AND 80;

CREATE OR REPLACE TEMP TABLE final_landmark_gates AS
SELECT *,
  ((deathtime IS NULL OR deathtime > landmark_time)
   AND (deathtime IS NOT NULL OR dod IS NULL
        OR dod >= CAST(landmark_time AS DATE))) AS alive_after_day5,
  outtime > landmark_time AS icu_after_day5,
  (antibiotic_time IS NOT NULL AND culture_time IS NOT NULL AND sofa_time IS NOT NULL
   AND antibiotic_time <= landmark_time AND culture_time <= landmark_time
   AND sofa_time <= landmark_time) AS sepsis_complete_by_day5
FROM final_sepsis_mother;

CREATE OR REPLACE TEMP TABLE final_magnesium_events AS
SELECT f.subject_id, f.hadm_id, f.stay_id,
       l.labevent_id, l.charttime, l.hadm_id AS source_hadm_id,
       l.valuenum AS raw_value, l.valueuom AS raw_unit,
       CASE WHEN lower(trim(l.valueuom)) = 'mg/dl' THEN CAST(l.valuenum AS DOUBLE)
            WHEN lower(trim(l.valueuom)) = 'mmol/l' THEN CAST(l.valuenum AS DOUBLE) * 2.4305
       END AS magnesium_mg_dl,
       CAST(floor(date_diff('microseconds', f.intime, l.charttime)
                  / 86400000000.0) + 1 AS INTEGER) AS day
FROM final_landmark_gates AS f
JOIN mimiciv_hosp.labevents AS l
  ON f.subject_id = l.subject_id
 AND (f.hadm_id = l.hadm_id OR l.hadm_id IS NULL)
WHERE l.itemid = 50960 AND l.valuenum > 0
  AND lower(trim(l.valueuom)) IN ('mg/dl', 'mmol/l')
  AND l.charttime >= f.intime AND l.charttime < f.landmark_time
  AND l.charttime < f.outtime;

CREATE OR REPLACE TEMP TABLE final_magnesium_daily AS
SELECT subject_id, hadm_id, stay_id, day,
       avg(magnesium_mg_dl) AS magnesium_mg_dl,
       count(*) AS raw_records
FROM final_magnesium_events
GROUP BY subject_id, hadm_id, stay_id, day;

CREATE OR REPLACE TEMP TABLE final_screened AS
WITH measurements AS (
  SELECT stay_id, count(*) AS n_observed_magnesium_days,
         sum(raw_records) AS n_magnesium_records,
         max(CASE WHEN day = 1 THEN magnesium_mg_dl END) AS mg_day1,
         max(CASE WHEN day = 2 THEN magnesium_mg_dl END) AS mg_day2,
         max(CASE WHEN day = 3 THEN magnesium_mg_dl END) AS mg_day3,
         max(CASE WHEN day = 4 THEN magnesium_mg_dl END) AS mg_day4,
         max(CASE WHEN day = 5 THEN magnesium_mg_dl END) AS mg_day5
  FROM final_magnesium_daily GROUP BY stay_id
)
SELECT f.*, coalesce(m.n_observed_magnesium_days, 0) AS n_observed_magnesium_days,
       coalesce(m.n_magnesium_records, 0) AS n_magnesium_records,
       m.mg_day1, m.mg_day2, m.mg_day3, m.mg_day4, m.mg_day5,
       CASE WHEN NOT f.alive_after_day5 THEN 'Death/date conflict by Day 5'
            WHEN NOT f.icu_after_day5 THEN 'ICU ended by Day 5'
            WHEN NOT f.sepsis_complete_by_day5 THEN 'Sepsis evidence incomplete by Day 5'
            WHEN coalesce(m.n_observed_magnesium_days, 0) < 3 THEN 'Fewer than 3 observed Mg windows'
            ELSE 'Included' END AS selection_reason
FROM final_landmark_gates AS f LEFT JOIN measurements AS m USING (stay_id);

CREATE OR REPLACE TEMP TABLE final_cohort AS
WITH deaths AS (
 SELECT *, coalesce(deathtime,
       CAST(dod AS TIMESTAMP) + INTERVAL '1' DAY - INTERVAL '1' MICROSECOND) AS death_time
 FROM final_screened WHERE selection_reason = 'Included'
), events AS (
 SELECT *, CAST(death_time IS NOT NULL AND death_time > landmark_time
                AND death_time <= endpoint_time AS INTEGER) AS event28
 FROM deaths
)
SELECT *, 5.0 AS entry_days,
       date_diff('microseconds', intime,
         CASE WHEN event28 = 1 THEN death_time ELSE endpoint_time END)
          / 86400000000.0 AS exit_days,
       CASE WHEN deathtime IS NOT NULL THEN 'Exact current-admission timestamp'
            WHEN dod IS NOT NULL THEN 'Date-only DOD, end-of-day rule'
            ELSE 'No death date recorded' END AS death_time_precision
FROM events;

CREATE OR REPLACE TEMP TABLE final_flow AS
SELECT 1 AS sequence, 'All ICU stays' AS stage, count(*) AS n FROM mimiciv_icu.icustays
UNION ALL SELECT 2, 'First ICU per patient', count(*) FROM final_first_icu
UNION ALL SELECT 3, 'First ICU and official Sepsis-3', count(*)
 FROM final_first_icu f JOIN mimiciv_derived.sepsis3 s USING(subject_id,stay_id) WHERE s.sepsis3 IS TRUE
UNION ALL SELECT 4, 'Also age 18-80 years', count(*) FROM final_screened
UNION ALL SELECT 5, 'Also alive after 120 h', count(*) FROM final_screened WHERE alive_after_day5
UNION ALL SELECT 6, 'Also ICU ongoing after 120 h', count(*) FROM final_screened WHERE alive_after_day5 AND icu_after_day5
UNION ALL SELECT 7, 'Also complete Sepsis-3 evidence by 120 h', count(*) FROM final_screened WHERE alive_after_day5 AND icu_after_day5 AND sepsis_complete_by_day5
UNION ALL SELECT 8, 'Also at least 3 observed Mg windows', count(*) FROM final_cohort;
