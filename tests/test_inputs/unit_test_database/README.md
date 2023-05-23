# Creating databases for unit tests

Create a local CAZyme database for the unit tests, containing the annotations for PL28 and CE5 CAZymes:

```bash
cazy_webscraper <email> -o tests/test_inputs/unit_test_database/unit_test_<DATE>.db --families PL28,CE5 -f -n
```

Then update the db path in `conftest.py`
