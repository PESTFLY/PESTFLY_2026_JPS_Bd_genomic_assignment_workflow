# Metadata and specimen identifiers

The public workbook `data/000_input_data/metadata.xlsx` has one sheet, `Selection`, with 352 rows and 13 columns in the supplied benchmark.

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

Step 01 requires `sample_id`, `is_reference`, `population` and `country`; downstream Step 02 also requires usable `macroregion_3` and `subregion` classifications. Do not silently replace them with geographic continent names. The benchmark Mascarene references use lineage based Asia and Southeast Asia coding, with alternative coding explicitly tested.

Step 01 cleans spaces, standardises column names and harmonises stated Congo and Reunion cases. It records missing labels and classifications rather than inferring all absent metadata. Allowed additional alignment labels are `Bdors` and `Blati`, the reference genome labels outside the specimen table.

Grouped geographic validation uses collection geography and site metadata to define held out units. The class being predicted remains the analytical reference class. Source collection and geography can be confounded; the audit reports representation without asserting that those effects have been removed statistically.

Declared commodity origin for intercepted larvae is separate reporting metadata. Step 07 supports an external table with `sample_id` and `sample_origin_country`, optional `sample_origin_country_source` and `commodity`. Its built in benchmark mapping is retained. Adult trap detections have no declared commodity country by default.
