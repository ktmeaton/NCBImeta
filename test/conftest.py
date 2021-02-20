"""
NCBImeta Conftest
This is the pytest fixture file.

@author: Katherine Eaton
"""

# -----------------------------------------------------------------------------#
#                         Modules and Packages                                 #
# -----------------------------------------------------------------------------#
import pytest

# -----------------------------------------------------------------------------#
#                             Fixtures                                         #
# -----------------------------------------------------------------------------#

# Reminder, awk was used to help efficiently format the column and value lists
# To get the column names:
# head -n 1 test_Assembly.txt | awk '
#     BEGIN{FS="\t"; ORS=",\n"}{for (i=1;i<=NF;i++){print "\x27"$i"\x27"}}'
# To get the metadata field values:
# tail -n+2 test_Assembly.txt | awk '
#    BEGIN{FS="\t"; ORS=", "}{for (i=1;i<=NF;i++){print "\x27"$i"\x27"}}'
# To format into single line lists, run through black pre-commit


@pytest.fixture(scope="module")
def assembly_table_data():
    """Return a dictionary containing the expected values of the Assembly Table"""

    columns = [
        "id",
        "Assembly_id",
        "AssemblyAccession",
        "AssemblyBioSampleAccession",
        "AssemblyGenbankBioprojectAccession",
        "AssemblyOrganism",
        "AssemblyContigCount",
        "AssemblyTotalLength",
        "AssemblySubmissionDate",
        "AssemblyFTPGenbank",
        "AssemblyComment",
    ]

    metadata = [
        "1",
        "5025191",
        "GCF_009295945.1",
        "SAMN12991206",
        "PRJNA269675",
        "Yersinia pestis (enterobacteria)",
        "5",
        "4759144",
        "2019/10/24 00:00",
        "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/295/945/GCA_009295945.1_ASM929594v1",
        "",
    ]

    table_dict = {}
    # Populate the dict with data
    for i in range(0, len(columns)):
        key = columns[i]
        value = metadata[i]
        table_dict[key] = value

    return table_dict


@pytest.fixture(scope="module")
def bioproject_table_data():
    """Return a dictionary containing the expected values of the BioProject Table"""

    columns = [
        "id",
        "BioProject_id",
        "BioProjectAccession",
        "BioProjectTitle",
        "BioProjectOrganismLabel",
        "BioProjectRegistrationDate",
        "BioProjectComment",
    ]

    metadata = [
        "1",
        "269675",
        "PRJNA269675",
        "Pathogenic microorganism Genome sequencing",
        "",
        "2014-12-09",
        "",
    ]

    table_dict = {}
    # Populate the dict with data
    for i in range(0, len(columns)):
        key = columns[i]
        value = metadata[i]
        table_dict[key] = value
    return table_dict


@pytest.fixture(scope="module")
def biosample_table_data():
    """Return a dictionary containing the expected values of the BioSample Table"""

    columns = [
        "id",
        "BioSample_id",
        "BioSampleAccession",
        "BioSampleAccessionSecondary",
        "BioSampleBioProjectAccession",
        "BioSampleSRAAccession",
        "BioSampleOrganism",
        "BioSampleStrain",
        "BioSampleSubmissionDate",
        "BioSampleComment",
    ]

    metadata = [
        "1",
        "12991206",
        "SAMN12991206",
        "",
        "",
        "SRS5502739",
        "TestOrganism1",
        "TestStrain1",
        "2019-10-08T07:15:03.953",
        "",
    ]

    table_dict = {}
    # Populate the dict with data
    for i in range(0, len(columns)):
        key = columns[i]
        value = metadata[i]
        table_dict[key] = value

    return table_dict


@pytest.fixture(scope="module")
def nucleotide_table_data():
    """Return a dictionary containing the expected values of the Nucleotide Table"""

    columns = [
        "id",
        "Nucleotide_id",
        "NucleotideAccession",
        "NucleotideBioProjectAccession",
        "NucleotideOrganism",
        "NucleotideLength",
        "NucleotideComment",
    ]

    metadata = [
        [
            "1",
            "1769169004",
            "NZ_CP045153",
            "PRJNA224116",
            "Yersinia pestis",
            "9610",
            "",
        ],
        [
            "2",
            "1769169003",
            "NZ_CP045152",
            "PRJNA224116",
            "Yersinia pestis",
            "33990",
            "",
        ],
        [
            "3",
            "1769169002",
            "NZ_CP045151",
            "PRJNA224116",
            "Yersinia pestis",
            "68343",
            "",
        ],
        [
            "4",
            "1769169001",
            "NZ_CP045150",
            "PRJNA224116",
            "Yersinia pestis",
            "100984",
            "",
        ],
        [
            "5",
            "1769169000",
            "NZ_CP045149",
            "PRJNA224116",
            "Yersinia pestis",
            "4546217",
            "",
        ],
        ["6", "1769094081", "CP045153", "PRJNA269675", "Yersinia pestis", "9610", ""],
        ["7", "1769094034", "CP045152", "PRJNA269675", "Yersinia pestis", "33990", ""],
        ["8", "1769093958", "CP045151", "PRJNA269675", "Yersinia pestis", "68343", ""],
        [
            "9",
            "1769093841",
            "CP045150",
            "PRJNA269675",
            "Yersinia pestis",
            "100984",
            "",
        ],
        [
            "10",
            "1769089848",
            "CP045149",
            "PRJNA269675",
            "Yersinia pestis",
            "4546217",
            "",
        ],
    ]

    table_dict = {}

    for c in range(0, len(columns)):
        col = columns[c]
        table_dict[col] = []
        for m in range(0, len(metadata)):
            print(col, c, metadata[m])
            meta = metadata[m][c]
            table_dict[col].append(meta)

    return table_dict


@pytest.fixture(scope="module")
def pubmed_table_data():
    """Return a dictionary containing the expected values of the Pubmed Table"""
    columns = [
        "id",
        "Pubmed_id",
        "PubmedPublishYear",
        "PubmedPubishMonth",
        "PubmedPublishDay",
        "PubmedEPublishDate",
        "PubmedPublishModel",
        "PubmedType",
        "PubmedComment",
    ]

    # Currently testing a pubmed file that is empty
    metadata = [
        "1",
        "26634751",
        "2015",
        "Dec",
        "03",
        "",
        "Electronic",
        "Journal Article",
        "",
    ]

    table_dict = {}
    # Populate the dict with data
    for i in range(0, len(columns)):
        key = columns[i]
        value = metadata[i]
        table_dict[key] = value

    return table_dict


@pytest.fixture(scope="module")
def sra_table_data():
    """Return a dictionary containing the expected values of the SRA Table"""

    columns = [
        "id",
        "SRA_id",
        "SRABioProjectAccession",
        "SRABioSampleAccession",
        "SRASampleAccession",
        "SRASampleName",
        "SRAExperimentAccession",
        "SRAExperimentName",
        "SRARunAccession",
        "SRAComment",
    ]

    metadata = [
        "1",
        "9179237",
        "PRJNA269675",
        "SAMN12991206",
        "SRS5502739",
        "2019 Sample 199",
        "SRX6977651",
        "2019 Sample 199 Library 02",
        "SRR10259780",
        "",
        "2",
        "9179236",
        "PRJNA269675",
        "SAMN12991206",
        "SRS5502739",
        "2019 Sample 199",
        "SRX6977650",
        "2019 Sample 199 Library 01",
        "SRR10259781",
        "",
    ]

    table_dict = {}
    # Populate the dict with data
    for metadata_i in range(0, len(metadata)):
        columns_i = metadata_i
        # If we're starting the next row of metadata
        if metadata_i >= len(columns):
            columns_i = metadata_i % len(columns)
        # Add key-value pair
        key = columns[columns_i]
        value = metadata[metadata_i]
        # Check if already present, if not add as list
        if key not in table_dict:
            table_dict[key] = [value]
        # Otherwise append, creating a list of values
        else:
            table_dict[key].append(value)

    return table_dict


@pytest.fixture(scope="module")
def assembly_xml():
    """Return a str of the expected xml from the assembly table"""
    assembly_xml_str = """
    <eSummaryResult>
        <DocumentSummarySet status="OK">
            <DbBuild>Build191210-0555.1</DbBuild>
            <DocumentSummary uid="5025191">
                <RsUid>14769648</RsUid>
                <GbUid>14768768</GbUid>
                <AssemblyAccession>GCF_009295945.1</AssemblyAccession>
                <LastMajorReleaseAccession>GCF_009295945.1</LastMajorReleaseAccession>
                <LatestAccession/>
                <ChainId>9295945</ChainId>
                <AssemblyName>ASM929594v1</AssemblyName>
                <UCSCName/>
                <EnsemblName/>
                <Taxid>632</Taxid>
                <Organism>Yersinia pestis (enterobacteria)</Organism>
                <SpeciesTaxid>632</SpeciesTaxid>
                <SpeciesName>Yersinia pestis</SpeciesName>
                <AssemblyType>haploid</AssemblyType>
                <AssemblyClass>haploid</AssemblyClass>
                <AssemblyStatus>Complete Genome</AssemblyStatus>
                <WGS/>
                <GB_BioProjects>
                    <Bioproj>
                        <BioprojectAccn>PRJNA269675</BioprojectAccn>
                        <BioprojectId>269675</BioprojectId>
                    </Bioproj>
                </GB_BioProjects>
                <GB_Projects>
                </GB_Projects>
                <RS_BioProjects>
                    <Bioproj>
                        <BioprojectAccn>PRJNA224116</BioprojectAccn>
                        <BioprojectId>224116</BioprojectId>
                    </Bioproj>
                </RS_BioProjects>
                <RS_Projects>
                </RS_Projects>
                <BioSampleAccn>SAMN12991206</BioSampleAccn>
                <BioSampleId>12991206</BioSampleId>
                <Biosource>
                    <InfraspeciesList>
                        <Infraspecie>
                            <Sub_type>strain</Sub_type>
                            <Sub_value>SCPM-O-DNA-18 (I-3113)</Sub_value>
                        </Infraspecie>
                    </InfraspeciesList>
                    <Sex/>
                    <Isolate/>
                </Biosource>
                <Coverage>176.1</Coverage>
                <PartialGenomeRepresentation>false</PartialGenomeRepresentation>
                <Primary>14769638</Primary>
                <AssemblyDescription/>
                <ReleaseLevel>Major</ReleaseLevel>
                <ReleaseType>Major</ReleaseType>
                <AsmReleaseDate_GenBank>2019/10/24 00:00</AsmReleaseDate_GenBank>
                <AsmReleaseDate_RefSeq>2019/10/25 00:00</AsmReleaseDate_RefSeq>
                <SeqReleaseDate>2019/10/24 00:00</SeqReleaseDate>
                <AsmUpdateDate>2019/10/25 00:00</AsmUpdateDate>
                <SubmissionDate>2019/10/24 00:00</SubmissionDate>
                <LastUpdateDate>2019/10/25 00:00</LastUpdateDate>
                <SubmitterOrganization>SRCAMB</SubmitterOrganization>
                <RefSeq_category>na</RefSeq_category>
                <AnomalousList>
                </AnomalousList>
                <ExclFromRefSeq>
                </ExclFromRefSeq>
                <PropertyList>
                    <string>full-genome-representation</string>
                    <string>genbank_has_annotation</string>
                    <string>has-chromosome</string>
                    <string>has-plasmid</string>
                    <string>has_annotation</string>
                    <string>latest</string>
                    <string>latest_genbank</string>
                    <string>latest_refseq</string>
                    <string>refseq_has_annotation</string>
                </PropertyList>
                <FromType/>
                <Synonym>
                    <Genbank>GCA_009295945.1</Genbank>
                    <RefSeq>GCF_009295945.1</RefSeq>
                    <Similarity>identical</Similarity>
                </Synonym>
                <ContigN50>4546217</ContigN50>
                <ScaffoldN50>4546217</ScaffoldN50>
                <FtpPath_GenBank>ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/295/945/GCA_009295945.1_ASM929594v1</FtpPath_GenBank>
                <FtpPath_RefSeq>ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/295/945/GCF_009295945.1_ASM929594v1</FtpPath_RefSeq>
                <FtpPath_Assembly_rpt>ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/295/945/GCF_009295945.1_ASM929594v1/GCF_009295945.1_ASM929594v1_assembly_report.txt</FtpPath_Assembly_rpt>
                <FtpPath_Stats_rpt>ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/295/945/GCF_009295945.1_ASM929594v1/GCF_009295945.1_ASM929594v1_assembly_stats.txt</FtpPath_Stats_rpt>
                <FtpPath_Regions_rpt/>
                <SortOrder>5C1X9954537820092959459898</SortOrder>
                <Meta>
                    <![CDATA[
                        <Stats>
                            <Stat category="alt_loci_count" sequence_tag="all">0</Stat>
                            <Stat category="chromosome_count" sequence_tag="all">1</Stat>
                            <Stat category="contig_count" sequence_tag="all">5</Stat>
                            <Stat category="contig_l50" sequence_tag="all">1</Stat>
                            <Stat category="contig_n50" sequence_tag="all">4546217</Stat>
                            <Stat category="non_chromosome_replicon_count" sequence_tag="all">4</Stat>
                            <Stat category="replicon_count" sequence_tag="all">5</Stat>
                            <Stat category="scaffold_count" sequence_tag="all">5</Stat>
                            <Stat category="scaffold_count" sequence_tag="placed">5</Stat>
                            <Stat category="scaffold_count" sequence_tag="unlocalized">0</Stat>
                            <Stat category="scaffold_count" sequence_tag="unplaced">0</Stat>
                            <Stat category="scaffold_l50" sequence_tag="all">1</Stat>
                            <Stat category="scaffold_n50" sequence_tag="all">4546217</Stat>
                            <Stat category="total_length" sequence_tag="all">4759144</Stat>
                            <Stat category="ungapped_length" sequence_tag="all">4759144</Stat>
                        </Stats>
                        <FtpSites>
                            <FtpPath type=Assembly_rpt>ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/295/945/GCF_009295945.1_ASM929594v1/GCF_009295945.1_ASM929594v1_assembly_report.txt</FtpPath>
                            <FtpPath type=GenBank>ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/295/945/GCA_009295945.1_ASM929594v1</FtpPath>
                            <FtpPath type=RefSeq>ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/295/945/GCF_009295945.1_ASM929594v1</FtpPath>
                            <FtpPath type=Stats_rpt>ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/295/945/GCF_009295945.1_ASM929594v1/GCF_009295945.1_ASM929594v1_assembly_stats.txt</FtpPath>
                        </FtpSites>
                        <assembly-level>5</assembly-level>
                        <assembly-status>Complete Genome</assembly-status>
                        <representative-status>na</representative-status>
                        <submitter-organization>SRCAMB</submitter-organization>
                    ]]>
                </Meta>
            </DocumentSummary>
        </DocumentSummarySet>
    </eSummaryResult>"""

    return assembly_xml_str
