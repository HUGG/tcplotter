from setuptools import setup

setup(name='tcplotter',
      version='0.2.1',
      description='Plots thermochronometer ages and closure temperatures',
      url='https://github.com/HUGG/tcplotter',
      author='David Whipp',
      author_email='david.whipp@helsinki.fi',
      license='MIT',
      packages=['tcplotter'],
      entry_points={
          'console_scripts': [
              'time-vs-temp=tcplotter.time_vs_temp:main',
              'eu-vs-radius=tcplotter.eu_vs_radius:main',
              'rate-vs-radius-eu=tcplotter.rate_vs_radius_eu:main',
              'rate-vs-age-tc=tcplotter.rate_vs_age_tc  :main',
          ],
      },
      zip_safe=False
      )
