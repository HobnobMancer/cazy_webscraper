##################################################################
# These unit tests are in the process of being updated
# This includes updating the paths to meet the new scraper structure
# This also includes factorising out the tests to make the script size
# easier to handle
#####################################################################


# #!/usr/bin/env python
# # -*- coding: utf-8 -*-
# # Author:
# # Emma E. M. Hobbs

# # Contact
# # eemh1@st-andrews.ac.uk

# # Emma E. M. Hobbs,
# # Biomolecular Sciences Building,
# # University of St Andrews,
# # North Haugh Campus,
# # St Andrews,
# # KY16 9ST
# # Scotland,
# # UK

# # The MIT License

# """Tests the module sql which builds and interacts with an SQL database.

# These test are intened to be run from the root of the repository using:
# pytest -v
# """


# import pytest

# from argparse import Namespace, ArgumentParser
# from pathlib import Path

# from scraper import expand
# from scraper.expand import get_pdb_structures, get_genbank_sequences
# from scraper.sql.sql_orm import Cazyme, Cazymes_Genbanks, Genbank, Taxonomy
# from scraper.utilities import file_io, parse_configuration


# @pytest.fixture
# def db_path():
#     path_ = Path("tests")
#     path_ = path_ / "test_inputs" / "test_inputs_expand" / "unit_test_2021-03-11--13-06-42.db"
#     return path_


# @pytest.fixture
# def output_dir(test_dir):
#     path_ = test_dir / "test_outputs"
#     return path_


# @pytest.fixture
# def tax_filter():
#     return set(["Nonlabens"])


# @pytest.fixture
# def genbank_query(db_session):
#     query = db_session.query(Genbank, Cazymes_Genbanks, Taxonomy).\
#         join(Genbank, (Genbank.genbank_id == Cazymes_Genbanks.genbank_id)).\
#         filter(Cazymes_Genbanks.primary == True).all()
#     return query


# # tests for get_pdb_structures


# def test_main_no_db(monkeypatch):
#     """Test main() when an the database file cannot be found."""

#     def mock_building_parser(*args, **kwargs):
#         parser_args = ArgumentParser(
#             prog="cazy_webscraper.py",
#             usage=None,
#             description="Scrape the CAZy database",
#             conflict_handler="error",
#             add_help=True,
#         )
#         return parser_args

#     def mock_parser(*args, **kwargs):
#         parser = Namespace(
#             database=Path("--"),
#             verbose=False,
#             log=None,
#             force=False,
#             nodelete=False,
#             outdir=None,
#         )
#         return parser

#     def mock_no_return(*args, **kwargs):
#         return

#     def mock_config(*args, **kwargs):
#         return None, set()

#     monkeypatch.setattr(utilities, "build_pdb_structures_parser", mock_building_parser)
#     monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
#     monkeypatch.setattr(parse_configuration, "get_configuration", mock_config)
#     monkeypatch.setattr(file_io, "make_output_directory", mock_no_return)

#     with pytest.raises(SystemExit) as pytest_wrapped_e:
#         get_pdb_structures.main()
#     assert pytest_wrapped_e.type == SystemExit


# def test_main_outdir_is_none(db_path, monkeypatch):
#     """Test main() when outdir=None."""

#     def mock_building_parser(*args, **kwargs):
#         parser_args = ArgumentParser(
#             prog="cazy_webscraper.py",
#             usage=None,
#             description="Scrape the CAZy database",
#             conflict_handler="error",
#             add_help=True,
#         )
#         return parser_args

#     def mock_parser(*args, **kwargs):
#         parser = Namespace(
#             database=db_path,
#             outdir=None,
#             verbose=False,
#             log=None,
#             force=False,
#             nodelete=False,
#         )
#         return parser

#     def mock_no_return(*args, **kwargs):
#         return

#     def mock_config(*args, **kwargs):
#         return None, set()

#     monkeypatch.setattr(utilities, "build_pdb_structures_parser", mock_building_parser)
#     monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
#     monkeypatch.setattr(utilities, "config_logger", mock_no_return)
#     monkeypatch.setattr(parse_configuration, "get_configuration", mock_config)
#     monkeypatch.setattr(file_io, "make_output_directory", mock_no_return)
#     monkeypatch.setattr(get_pdb_structures, "get_every_cazymes_structures", mock_no_return)

#     get_pdb_structures.main()


# def test_main(db_path, output_dir, monkeypatch):
#     """Test main()."""

#     def mock_building_parser(*args, **kwargs):
#         parser_args = ArgumentParser(
#             prog="cazy_webscraper.py",
#             usage=None,
#             description="Scrape the CAZy database",
#             conflict_handler="error",
#             add_help=True,
#         )
#         return parser_args

#     def mock_parser(*args, **kwargs):
#         parser = Namespace(
#             database=db_path,
#             outdir=output_dir,
#             verbose=False,
#             log=None,
#             force=True,
#             nodelete=True,
#         )
#         return parser

#     def mock_no_return(*args, **kwargs):
#         return

#     def mock_config(*args, **kwargs):
#         return None, set()

#     monkeypatch.setattr(utilities, "build_pdb_structures_parser", mock_building_parser)
#     monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
#     monkeypatch.setattr(utilities, "config_logger", mock_no_return)
#     monkeypatch.setattr(parse_configuration, "get_configuration", mock_config)
#     monkeypatch.setattr(file_io, "make_output_directory", mock_no_return)
#     monkeypatch.setattr(get_pdb_structures, "get_every_cazymes_structures", mock_no_return)

#     get_pdb_structures.main()


# def test_main_argv(db_path, output_dir, monkeypatch):
#     """Test main()."""

#     def mock_building_parser(*args, **kwargs):
#         parser_args = ArgumentParser(
#             prog="cazy_webscraper.py",
#             usage=None,
#             description="Scrape the CAZy database",
#             conflict_handler="error",
#             add_help=True,
#         )
#         return parser_args

#     def mock_parser(*args, **kwargs):
#         parser = Namespace(
#             database=db_path,
#             outdir=output_dir,
#             verbose=False,
#             log=None,
#             force=True,
#             nodelete=True,
#         )
#         return parser

#     def mock_no_return(*args, **kwargs):
#         return

#     def mock_config(*args, **kwargs):
#         return {}, set()

#     monkeypatch.setattr(utilities, "build_pdb_structures_parser", mock_building_parser)
#     monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
#     monkeypatch.setattr(utilities, "config_logger", mock_no_return)
#     monkeypatch.setattr(parse_configuration, "get_configuration", mock_config)
#     monkeypatch.setattr(file_io, "make_output_directory", mock_no_return)
#     monkeypatch.setattr(get_pdb_structures, "get_structures_for_specific_cazymes", mock_no_return)

#     get_pdb_structures.main()


# def test_get_every_cazymes_structures_primary(db_session, output_dir, monkeypatch):
#     """Test get_every_cazymes_structures() when primary is True and taxonomy filter is given."""

#     def mock_no_return(*args, **kwargs):
#         return

#     monkeypatch.setattr(get_pdb_structures, "download_pdb_structures", mock_no_return)

#     args = {"args": Namespace(primary=True)}
#     tax_filter = set(["Nonlabens"])

#     get_pdb_structures.get_every_cazymes_structures(
#         output_dir,
#         tax_filter,
#         db_session,
#         args["args"],
#     )


# def test_get_every_cazymes_structures_all(db_session, output_dir, monkeypatch):
#     """Test get_every_cazymes_structures() when primary is False and no taxonomy filter is given."""

#     def mock_no_return(*args, **kwargs):
#         return

#     monkeypatch.setattr(get_pdb_structures, "download_pdb_structures", mock_no_return)

#     args = {"args": Namespace(primary=False)}
#     tax_filter = None

#     get_pdb_structures.get_every_cazymes_structures(
#         output_dir,
#         tax_filter,
#         db_session,
#         args["args"],
#     )


# def test_get_structures_for_specific_cazymes_primary(db_session, output_dir, monkeypatch):
#     """Test get_structures_for_specific_cazymes when primary is true and tax filter is given."""

#     def mock_no_return(*args, **kwargs):
#         return

#     monkeypatch.setattr(get_pdb_structures, "download_pdb_structures", mock_no_return)

#     args = {"args": Namespace(primary=True)}
#     tax_filter = set(["Nonlabens"])
#     config_dict = {"classes": ["PL"], "Polysaccharide Lyases (PLs)": ["PL28"]}

#     get_pdb_structures.get_structures_for_specific_cazymes(
#         output_dir,
#         config_dict,
#         tax_filter,
#         db_session,
#         args["args"],
#     )


# def test_get_structures_for_specific_cazymes(db_session, output_dir, monkeypatch):
#     """Test get_structures_for_specific_cazymes when primary is False and tax filter not given."""

#     def mock_no_return(*args, **kwargs):
#         return

#     monkeypatch.setattr(get_pdb_structures, "download_pdb_structures", mock_no_return)

#     args = {"args": Namespace(primary=False)}
#     tax_filter = None
#     config_dict = {
#         "classes": ["PL"],
#         "Polysaccharide Lyases (PLs)": ["PL28","GH3_1"],
#         "CAZyclass": None,
#     }

#     get_pdb_structures.get_structures_for_specific_cazymes(
#         output_dir,
#         config_dict,
#         tax_filter,
#         db_session,
#         args["args"],
#     )


# # test for get_genbank_sequences


# def test_main_no_db_genbank(monkeypatch):
#     """Test main() when an the database file cannot be found."""

#     def mock_building_parser(*args, **kwargs):
#         parser_args = ArgumentParser(
#             prog="cazy_webscraper.py",
#             usage=None,
#             description="Scrape the CAZy database",
#             conflict_handler="error",
#             add_help=True,
#         )
#         return parser_args

#     def mock_parser(*args, **kwargs):
#         parser = Namespace(
#             database=Path("--"),
#             email="dummy_email",
#             verbose=False,
#             log=None,
#             force=False,
#             nodelete=False,
#             outdir=None,
#         )
#         return parser

#     def mock_no_return(*args, **kwargs):
#         return

#     def mock_config(*args, **kwargs):
#         return None, set()

#     monkeypatch.setattr(utilities, "build_genbank_sequences_parser", mock_building_parser)
#     monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
#     monkeypatch.setattr(parse_configuration, "get_configuration", mock_config)
#     monkeypatch.setattr(file_io, "make_output_directory", mock_no_return)

#     with pytest.raises(SystemExit) as pytest_wrapped_e:
#         get_genbank_sequences.main()
#     assert pytest_wrapped_e.type == SystemExit


# def test_main_no_config_no_update(db_path, monkeypatch):
#     """Test main() when outdir=None."""

#     def mock_building_parser(*args, **kwargs):
#         parser_args = ArgumentParser(
#             prog="cazy_webscraper.py",
#             usage=None,
#             description="Scrape the CAZy database",
#             conflict_handler="error",
#             add_help=True,
#         )
#         return parser_args

#     def mock_parser(*args, **kwargs):
#         parser = Namespace(
#             database=db_path,
#             email="dummy_email",
#             outdir=None,
#             verbose=False,
#             log=None,
#             force=False,
#             nodelete=False,
#             update=False,
#         )
#         return parser

#     def mock_no_return(*args, **kwargs):
#         return

#     def mock_config(*args, **kwargs):
#         return None, set()

#     monkeypatch.setattr(utilities, "build_genbank_sequences_parser", mock_building_parser)
#     monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
#     monkeypatch.setattr(utilities, "config_logger", mock_no_return)
#     monkeypatch.setattr(parse_configuration, "get_configuration", mock_config)
#     monkeypatch.setattr(file_io, "make_output_directory", mock_no_return)
#     monkeypatch.setattr(
#         get_genbank_sequences,
#         "get_missing_sequences_for_everything",
#         mock_no_return,
#     )

#     get_genbank_sequences.main()


# def test_main_no_config_yes_update(db_path, monkeypatch):
#     """Test main() when outdir=None."""

#     def mock_building_parser(*args, **kwargs):
#         parser_args = ArgumentParser(
#             prog="cazy_webscraper.py",
#             usage=None,
#             description="Scrape the CAZy database",
#             conflict_handler="error",
#             add_help=True,
#         )
#         return parser_args

#     def mock_parser(*args, **kwargs):
#         parser = Namespace(
#             database=db_path,
#             email="dummy_email",
#             outdir=None,
#             verbose=False,
#             log=None,
#             force=False,
#             nodelete=False,
#             update=True,
#         )
#         return parser

#     def mock_no_return(*args, **kwargs):
#         return

#     def mock_config(*args, **kwargs):
#         return None, set()

#     monkeypatch.setattr(utilities, "build_genbank_sequences_parser", mock_building_parser)
#     monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
#     monkeypatch.setattr(utilities, "config_logger", mock_no_return)
#     monkeypatch.setattr(parse_configuration, "get_configuration", mock_config)
#     monkeypatch.setattr(file_io, "make_output_directory", mock_no_return)
#     monkeypatch.setattr(
#         get_genbank_sequences,
#         "add_and_update_all_sequences",
#         mock_no_return,
#     )

#     get_genbank_sequences.main()


# def test_main_yes_config_no_update(db_path, monkeypatch):
#     """Test main() when outdir=None."""

#     def mock_building_parser(*args, **kwargs):
#         parser_args = ArgumentParser(
#             prog="cazy_webscraper.py",
#             usage=None,
#             description="Scrape the CAZy database",
#             conflict_handler="error",
#             add_help=True,
#         )
#         return parser_args

#     def mock_parser(*args, **kwargs):
#         parser = Namespace(
#             database=db_path,
#             email="dummy_email",
#             outdir=None,
#             verbose=False,
#             log=None,
#             force=False,
#             nodelete=False,
#             update=False,
#         )
#         return parser

#     def mock_no_return(*args, **kwargs):
#         return

#     def mock_config(*args, **kwargs):
#         return {}, set()

#     monkeypatch.setattr(utilities, "build_genbank_sequences_parser", mock_building_parser)
#     monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
#     monkeypatch.setattr(utilities, "config_logger", mock_no_return)
#     monkeypatch.setattr(parse_configuration, "get_configuration", mock_config)
#     monkeypatch.setattr(file_io, "make_output_directory", mock_no_return)
#     monkeypatch.setattr(
#         get_genbank_sequences,
#         "get_missing_sequences_for_specific_records",
#         mock_no_return,
#     )

#     get_genbank_sequences.main()


# def test_main_yes_config_yes_update(db_path, monkeypatch):
#     """Test main() when outdir=None."""

#     def mock_building_parser(*args, **kwargs):
#         parser_args = ArgumentParser(
#             prog="cazy_webscraper.py",
#             usage=None,
#             description="Scrape the CAZy database",
#             conflict_handler="error",
#             add_help=True,
#         )
#         return parser_args

#     def mock_parser(*args, **kwargs):
#         parser = Namespace(
#             database=db_path,
#             email="dummy_email",
#             outdir=None,
#             verbose=False,
#             log=None,
#             force=False,
#             nodelete=False,
#             update=True,
#         )
#         return parser

#     def mock_no_return(*args, **kwargs):
#         return

#     def mock_config(*args, **kwargs):
#         return {}, set()

#     monkeypatch.setattr(utilities, "build_genbank_sequences_parser", mock_building_parser)
#     monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
#     monkeypatch.setattr(utilities, "config_logger", mock_no_return)
#     monkeypatch.setattr(parse_configuration, "get_configuration", mock_config)
#     monkeypatch.setattr(file_io, "make_output_directory", mock_no_return)
#     monkeypatch.setattr(
#         get_genbank_sequences,
#         "update_sequences_for_specific_records",
#         mock_no_return,
#     )

#     get_genbank_sequences.main()


# def test_get_missing_sequences_for_everything_primary(db_session, monkeypatch):
#     """Tests get_missing_sequences_for_everything() when primary is True and tax filter is given."""

#     def mock_accession(*args, **kwargs):
#         return []

#     def mock_no_return(*args, **kwargs):
#         return

#     monkeypatch.setattr(get_genbank_sequences, "extract_accessions", mock_accession)
#     monkeypatch.setattr(get_genbank_sequences, "get_sequences_add_to_db", mock_no_return)

#     args = {"args": Namespace(primary=True)}
#     tax_filter = set(["Nonlabens"])

#     get_genbank_sequences.get_missing_sequences_for_everything(
#         "date",
#         tax_filter,
#         db_session,
#         args["args"],
#     )


# def test_get_missing_sequences_for_everything(db_session, monkeypatch):
#     """Tests get_missing_sequences_for_everything() when primary is False and tax filter is None."""

#     def mock_accession(*args, **kwargs):
#         return ["acc1", "acc2"]

#     def mock_no_return(*args, **kwargs):
#         return

#     monkeypatch.setattr(get_genbank_sequences, "extract_accessions", mock_accession)
#     monkeypatch.setattr(get_genbank_sequences, "get_accession_chunks", mock_accession)
#     monkeypatch.setattr(get_genbank_sequences, "get_sequences_add_to_db", mock_no_return)

#     args = {"args": Namespace(primary=True, epost=200)}
#     tax_filter = None

#     get_genbank_sequences.get_missing_sequences_for_everything(
#         "date",
#         tax_filter,
#         db_session,
#         args["args"],
#     )


# def test_add_and_update_all_sequences_primary(db_session, monkeypatch):
#     """Tests add_and_update_all_sequences() when primary is True and tax filter is given."""

#     def mock_accession(*args, **kwargs):
#         return {}

#     monkeypatch.setattr(get_genbank_sequences, "extract_accessions_and_dates", mock_accession)

#     args = {"args": Namespace(primary=True, epost=200)}
#     tax_filter = set(["Nonlabens"])

#     get_genbank_sequences.add_and_update_all_sequences(
#         "date",
#         tax_filter,
#         db_session,
#         args["args"],
#     )


# def test_add_and_update_all_sequences_no_updates(db_session, monkeypatch):
#     """Tests get_missing_sequences_for_everything() when there are no seq to update."""

#     def mock_accession(*args, **kwargs):
#         return {"acc1":0, "acc2":1}

#     def mock_no_acc(*args, **kwargs):
#         return []

#     monkeypatch.setattr(get_genbank_sequences, "extract_accessions_and_dates", mock_accession)
#     monkeypatch.setattr(get_genbank_sequences, "get_accessions_for_new_sequences", mock_no_acc)

#     args = {"args": Namespace(primary=True, epost=200)}
#     tax_filter = None

#     get_genbank_sequences.add_and_update_all_sequences(
#         "date",
#         tax_filter,
#         db_session,
#         args["args"],
#     )


# def test_add_and_update_all_sequences(db_session, monkeypatch):
#     """Tests get_missing_sequences_for_everything() when primary is False and tax filter is None."""

#     def mock_accession(*args, **kwargs):
#         return {"acc1":1, "acc2":2}

#     def mock_acc(*args, **kwargs):
#         return ["acc", "acc1"]

#     def mock_no_return(*args, **kwargs):
#         return

#     monkeypatch.setattr(get_genbank_sequences, "extract_accessions_and_dates", mock_accession)
#     monkeypatch.setattr(get_genbank_sequences, "get_accessions_for_new_sequences", mock_acc)
#     monkeypatch.setattr(get_genbank_sequences, "get_accession_chunks", mock_acc)
#     monkeypatch.setattr(get_genbank_sequences, "get_sequences_add_to_db", mock_no_return)

#     args = {"args": Namespace(primary=True, epost=200)}
#     tax_filter = None

#     get_genbank_sequences.add_and_update_all_sequences(
#         "date",
#         tax_filter,
#         db_session,
#         args["args"],
#     )


# def test_get_missing_sequences_for_specific_records_primary(db_session, monkeypatch):
#     """Tests get_missing_sequences_for_specific_records, primary is True, tax filter isn't None."""

#     def mock_accession(*args, **kwargs):
#         return []

#     def mock_no_return(*args, **kwargs):
#         return

#     monkeypatch.setattr(get_genbank_sequences, "extract_accessions", mock_accession)
#     monkeypatch.setattr(get_genbank_sequences, "get_sequences_add_to_db", mock_no_return)

#     args = {"args": Namespace(primary=True, epost=150)}
#     tax_filter = None
#     config = {"classes":["PL"], "PL": ["PL28", "PL29_1"], "GH": None}

#     get_genbank_sequences.get_missing_sequences_for_specific_records(
#         "date",
#         config,
#         tax_filter,
#         db_session,
#         args["args"],
#     )


# def test_get_missing_sequences_for_specific_records(db_session, monkeypatch):
#     """Tests get_missing_sequences_for_specific_records, primary is False, tax filter isn't None."""

#     def mock_accession(*args, **kwargs):
#         return ["acc", "acc1"]

#     def mock_no_return(*args, **kwargs):
#         return

#     monkeypatch.setattr(get_genbank_sequences, "extract_accessions", mock_accession)
#     monkeypatch.setattr(get_genbank_sequences, "get_sequences_add_to_db", mock_no_return)

#     args = {"args": Namespace(primary=False, epost=150)}
#     tax_filter = set(["Nonlabens"])
#     config = {"classes":["PL"], "PL": ["PL28", "PL29_1"], "GH": None}

#     get_genbank_sequences.get_missing_sequences_for_specific_records(
#         "date",
#         config,
#         tax_filter,
#         db_session,
#         args["args"],
#     )


# def test_update_sequences_for_specific_records_primary_no_acc1(db_session, monkeypatch):
#     """Test update_sequences_for_specific_records, primary is True, tax filter not given."""

#     def mock_accession(*args, **kwargs):
#         return {}

#     def mock_acc(*args, **kwargs):
#         return []

#     def mock_no_return(*args, **kwargs):
#         return

#     monkeypatch.setattr(get_genbank_sequences, "extract_accessions_and_dates", mock_accession)
#     monkeypatch.setattr(get_genbank_sequences, "get_accessions_for_new_sequences", mock_acc)
#     monkeypatch.setattr(get_genbank_sequences, "get_accession_chunks", mock_acc)
#     monkeypatch.setattr(get_genbank_sequences, "get_sequences_add_to_db", mock_no_return)

#     args = {"args": Namespace(primary=True, epost=200)}
#     tax_filter = None
#     config = {"classes": ["PL"], "PL": ["PL28", "PL29_1"], "GH": None}

#     get_genbank_sequences.update_sequences_for_specific_records(
#         "date",
#         config,
#         tax_filter,
#         db_session,
#         args["args"],
#     )


# def test_update_sequences_for_specific_records_primary_no_acc2(db_session, monkeypatch):
#     """Test update_sequences_for_specific_records, primary is True, tax filter not given."""

#     def mock_accession(*args, **kwargs):
#         return {"acc1": 1, "acc2": 2}

#     def mock_acc(*args, **kwargs):
#         return []

#     def mock_no_return(*args, **kwargs):
#         return

#     monkeypatch.setattr(get_genbank_sequences, "extract_accessions_and_dates", mock_accession)
#     monkeypatch.setattr(get_genbank_sequences, "get_accessions_for_new_sequences", mock_acc)
#     monkeypatch.setattr(get_genbank_sequences, "get_accession_chunks", mock_acc)
#     monkeypatch.setattr(get_genbank_sequences, "get_sequences_add_to_db", mock_no_return)

#     args = {"args": Namespace(primary=True, epost=200)}
#     tax_filter = None
#     config = {"classes": ["PL"], "PL": ["PL28", "PL29_1"], "GH": None}

#     get_genbank_sequences.update_sequences_for_specific_records(
#         "date",
#         config,
#         tax_filter,
#         db_session,
#         args["args"],
#     )


# def test_update_sequences_for_specific_records_primary_false(db_session, monkeypatch):
#     """Test update_sequences_for_specific_records, primary is False, tax filter given."""

#     def mock_accession(*args, **kwargs):
#         return {"acc1": 1, "acc2": 2}

#     def mock_acc(*args, **kwargs):
#         return ["acc", "acc1"]

#     def mock_no_return(*args, **kwargs):
#         return

#     monkeypatch.setattr(get_genbank_sequences, "extract_accessions_and_dates", mock_accession)
#     monkeypatch.setattr(get_genbank_sequences, "get_accessions_for_new_sequences", mock_acc)
#     monkeypatch.setattr(get_genbank_sequences, "get_accession_chunks", mock_acc)
#     monkeypatch.setattr(get_genbank_sequences, "get_sequences_add_to_db", mock_no_return)

#     args = {"args": Namespace(primary=False, epost=200)}
#     tax_filter = set(["Nonlabens"])
#     config = {"classes": ["PL"], "PL": ["PL28", "PL29_1"], "GH": None}

#     get_genbank_sequences.update_sequences_for_specific_records(
#         "date",
#         config,
#         tax_filter,
#         db_session,
#         args["args"],
#     )


# def test_extract_accessions_no_tax(genbank_query):
#     """Test extract_accessions() when no tax filter is given."""
#     get_genbank_sequences.extract_accessions(genbank_query, None)


# def test_extract_accessions_tax_given(genbank_query, tax_filter):
#     """Test extract_accessions() when no tax filter is given."""
#     get_genbank_sequences.extract_accessions(genbank_query, tax_filter)


# def test_extract_accessions_and_dates_no_tax(genbank_query):
#     """Test extract_accessions() when no tax filter is given."""
#     get_genbank_sequences.extract_accessions_and_dates(genbank_query, None)


# def test_extract_accessions_and_dates_tax_given(genbank_query, tax_filter):
#     """Test extract_accessions() when no tax filter is given."""
#     get_genbank_sequences.extract_accessions_and_dates(genbank_query, tax_filter)


# def test_get_accession_chunks():
#     """Test get_accession_chunks()"""
#     lst = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
#     get_genbank_sequences.get_accession_chunks(lst, 2)


# def test_entry_retry():
#     """Test entrez_retry."""

#     def mock_record(*args, **kwargs):
#         return "test_record"

#     assert "test_record" == get_genbank_sequences.entrez_retry(mock_record)


# def test_entrez_retry_none():
#     """Test entrez_retry when nothing is returned."""

#     def mock_record(*args, **kwargs):
#         return

#     assert get_genbank_sequences.entrez_retry(mock_record) is None 
