
from setuptools import setup, find_namespace_packages

setup(
    name='bifrost_reporter',
    version='0.1.0',
    packages=find_namespace_packages(include=['bifrost_reporter',
                                              'bifrost_reporter.*'
                                             ]),
    install_requires=['pandas',
                      'pymongo',
                      'envyaml'],
    include_package_data=True,
    zip_safe=False,
    author='Your Name',
    author_email='your.email@example.com',
    description='A description of your package',
    url='https://github.com/yourusername/bifrost_reporter',
)