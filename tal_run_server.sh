export PYTHONPATH=`pwd`:`pwd`/charge_backend

FLASK_INCHI_DB=azf/inchi_mapping.json LIVAI_BASE_URL="https://livai-api.llnl.gov/v1" LIVAI_API_KEY="sk-KmXE6B6GqWz4IH6mLkGHog" FLASK_STOCK_DB=azf/zinc_stock.hdf5 python charge_backend/charge_server.py --port 8001 --host 0.0.0.0 --model gpt-5.1 --backend livai --config-file azf/config.yml

#(.venv) ) steady10:flask-copilot learsch1$ FLASK_INCHI_DB=azf/inchi_mapping.json LIVASE_URL="https://livai-api.llnl.gov/v1" LIVAI_API_KEY="sk-KmXE6B6GqWz4IH6mLkGHog" FLASK_STOCK_DB=azf/zinc_stock.hdf5 python charge_backend/charge_server.py --port 8001 --host 0.0.0.0 --model gpt-5.1 --backend livai --config-file azf/config.yml