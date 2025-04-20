from setuptools import setup, Extension

setup(
    name="mykmeanssp",
    version="0.0.4",
    author="Maha Massalha",
    author_email="mass.maha19@gmail.com",
    license="BSD",
    ext_modules=[
        Extension(
            'mykmeanssp',
            ['kmeans.c'],
        )
    ]
)
