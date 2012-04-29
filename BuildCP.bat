@echo off

python -c "import shutil; shutil.rmtree('build',ignore_errors=True)"
python setup.py --verbose build
REM ~ python setup.py bdist_msi
python setup.py install
python setup.py clean

