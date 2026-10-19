-- Stage 0 only: inspect public NCBI Pathogen Detection table schemas and
-- estimate the available meropenem cohort. Run in BigQuery SQL Workspace.
-- NCBI documents these tables as ncbi-pathogen-detect.pdbrowser.isolates
-- and .ast. The AST feed is an alpha release and updates daily.

-- 1) Inspect fields before relying on any field names in later queries.
SELECT table_name, column_name, data_type
FROM `ncbi-pathogen-detect.pdbrowser.INFORMATION_SCHEMA.COLUMNS`
WHERE table_name IN ('isolates', 'ast')
ORDER BY table_name, ordinal_position;

-- 2) NCBI-documented style of antibiotic counts, restricted to the target
-- species. Check the returned taxgroup_name spelling before proceeding.
SELECT taxgroup_name, COUNT(DISTINCT target_acc) AS isolates_tested
FROM `ncbi-pathogen-detect.pdbrowser.ast`
WHERE LOWER(antibiotic) = 'meropenem'
  AND LOWER(taxgroup_name) LIKE '%klebsiella pneumoniae%'
GROUP BY taxgroup_name
ORDER BY isolates_tested DESC;

-- 3) Once the schema is confirmed, use this join to count unique isolates
-- with meropenem AST and isolate-browser genotype records. Inspect source
-- type and AST interpretation fields before deciding the clinical subset
-- and resistance categories. Do not export data until the cohort is feasible.
SELECT
  COUNT(DISTINCT a.target_acc) AS matched_isolates,
  COUNT(DISTINCT IF(i.amr_genotypes IS NOT NULL
                    AND i.amr_genotypes != '', a.target_acc, NULL)) AS with_genotype_text
FROM `ncbi-pathogen-detect.pdbrowser.ast` AS a
JOIN `ncbi-pathogen-detect.pdbrowser.isolates` AS i
  USING (target_acc)
WHERE LOWER(a.antibiotic) = 'meropenem'
  AND LOWER(a.taxgroup_name) LIKE '%klebsiella pneumoniae%';

-- Note: if the schema uses a different genotype column name, update the
-- reference above after reviewing query 1. To define a clinical cohort,
-- inspect the isolate table's source/isolation fields and their values;
-- do not assume one field or its labels without checking.
