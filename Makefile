.PHONY : clean build upload

clean:
	rm -rf build
	rm -rf dist
	rm -rf cellatlas.egg-info
	rm -rf docs/_build
	rm -rf docs/api
	rm -rf .coverage

build:
	python -m build --wheel

upload:
	twine upload dist/*