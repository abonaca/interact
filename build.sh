rm -r build
rm *so*

python3 setup.py build_ext --inplace
