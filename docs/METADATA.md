# Metadata and specimen identifiers

`data/000_input_data/metadata.xlsx`, sheet `Selection`, contains 352 specimens and 13 columns.

| Original column | Role |
| :--- | :--- |
| `sample_id` | Unique specimen identifier matching alignment labels |
| `is_reference` | Reference TRUE or query FALSE; 330 TRUE and 22 FALSE |
| `macroregion` | Retained source metadata field; downstream candidate classes use `macroregion_3` |
| `population` | Population grouping used during QC and auditing |
| `country` | Recorded collection country; separate from declared commodity origin |
| `Site` | Collection site; names are standardised to lower case by Step 01 |
| `collection` | Source collection or study annotation |
| `sample_accession` | Associated accession retained for provenance and auditing |
| `organism` | Recorded taxonomic identification |
| `collection_date` | Recorded specimen collection date |
| `geo_loc_name` | Additional locality metadata |
| `macroregion_3` | Analytical macroregion class: Africa, Asia or Other in the benchmark |
| `subregion` | Conditional analytical subregion class |

Step 01 requires `sample_id`, `is_reference`, `population` and `country`. Step 02 also needs `macroregion_3` and `subregion`. Retain the supplied analytical classifications when reproducing the benchmark; Mascarene references use Asia and Southeast Asia coding.

Step 01 standardises column names and spaces, harmonises the stated Congo and Reunion cases, and records label and metadata issues. `Bdors` and `Blati` are allowed reference genome alignment labels outside the specimen table.

Grouped validation uses collection country and site to define holdouts while predicting the analytical class. Declared commodity origin is separate metadata: Step 07 accepts `sample_id` and `sample_origin_country`, with optional `sample_origin_country_source` and `commodity`. It includes the benchmark larval mapping; adult trap detections have no declared commodity country by default.
