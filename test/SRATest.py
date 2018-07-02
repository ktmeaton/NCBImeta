OUTPUT_DIR = "../NCBImeta_ouput/"
EMAIL = "ktmeaton@gmail.com"
DATABASE = "yersinia_pestis_db.sqlite"
TABLES = ["SRA"]
SEARCH_TERMS = {"SRA": "Yersinia pestis[Orgn]"}

TABLE_COLUMNS = {"SRA": [
                    {"CreateDate" : "CreateDate"},
                    {"UpdateDate" : "UpdateDate"},
                    {"Title" : "Title"},
                    {"Platform" : "Platform"},
                    {"LibraryName" : "LIBRARY_NAME"},
                    {"LibraryStrategy" : "LIBRARY_STRATEGY"},
                    {"LibrarySource" : "LIBRARY_SOURCE"},
                    {"LibrarySelection" : "LIBRARY_SELECTION"},
                    {"LibraryLayout" : "LIBRARY_LAYOUT"},
                    {"BioProjectAccession" : "Bioproject"},
                    {"BioSampleAccession" : "Biosample"},
                    {"ClusterName" : ["Statistics", "cluster_name"]},
                    {"TotalBases" : ["Statistics", "total_bases"]},
                    {"TotalSpots" : ["Statistics", "total_spots"]},
                    {"TotalRuns" : ["Statistics", "total_runs"]},
                    {"TotalSize" : ["Statistics", "total_size"]},
                    {"StaticDataAvailable" : ["Run", "static_data_available"]},
                    {"IsPublic" : ["Run", "is_public"]},
                    {"InstrumentModel" : ["Platform", "instrument_model"]},
                    {"SubmitterAccession" : ["Submitter", "acc"]},
                    {"CenterName" : ["Submitter", "center_name"]},
                    {"ContactName" : ["Submitter", "contact_name"]},
                    {"LabName" : ["Submitter", "lab_name"]},
                    {"ExperimentAccession" : ["Experiment", "acc"]},
                    {"ExperimentName" : ["Experiment", "name"]},
                    {"ExperimentStatus" : ["Experiment", "status"]},
                    {"ExperimentVersion" : ["Experiment", "ver"]},
                    {"StudyAcc" : ["Study", "acc"]},
                    {"StudName" : ["Study", "name"]},
                    {"OrganismName" : ["Organism", "CommonName"]},
                    {"OrganismTaxID" : ["Organism", "taxid"]},
                    {"SampleAccession" : ["Sample", "acc"]},
                    {"SampleName" : ["Sample", "name"]},
                    {"RunAccession" : ["Run", "acc"]}
                    ]}