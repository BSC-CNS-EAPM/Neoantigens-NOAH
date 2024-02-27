import setuptools

setuptools.setup(
    name='noah',
    version='0.1',
    packages=setuptools.find_packages(),
    url='https://github.com/BSC-CNS-EAPM/Neoantigens-NOAH',
    license='MIT',
    author='Pep Amengual, Oriol Gracia, Roc Farriol, Albert Cañellas',
    author_email='',
    description='A python package of NOAH',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    install_requires=[
        # list of packages your project depends on
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
    ],
    python_requires='>=3.6',
)