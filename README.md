# swissrepeat
Collaboration with ZHAW on annotation and analysis of the SwissProt database in terms of tandem repeats and disorder.

The code is split in 1+3 parts.
It's started with the data collection, analysis and processing.
The other three parts are for the evaluation and exploration grouped into three Scripts:
swissrepeats.Rmd: which analyses the repeat data (analogous to marcotte et al.)
swissdisorder.Rmd: foo (TODO)
empirical_and_expected_homorepeat_counts.Rmd: foo (TODO)
ordered_and_disorderd_swissrepeats.Rmd: foo (TODO)
and finally:
manuscript_figures_and_data.Rmd: which contains the figures and data from the paper.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them

```
Give examples
```

### Installing
Make a copy of the ```local_config_TEMPLATE.R``` file by:
```
cp local_config_TEMPLATE.R local_config.R
```
In the newly created ```local_config.R``` specify system specific configurations.

Install the required R-packages:
```
install.packages(c("plyr", "pals", "ggplot2", "scales", "grid", "RColorBrewer", "tidyverse"))
```

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
