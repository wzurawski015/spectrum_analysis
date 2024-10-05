from setuptools import setup

setup(
    name='spectrum_analysis_tool',
    version='1.02',
    description='Narzędzie do analizy spektralnej danych autokorelacyjnych',
    author='Twoje Imię',
    packages=['spectrum_analysis'],
    install_requires=[
        'numpy',
        'matplotlib',
        'pyfftw',
        'plotly',
        'yaml',
        'scipy'
    ],
    entry_points={
        'console_scripts': [
            'spectrum_analysis=spectrum_analysis_optimized:main'
        ]
    },
)
