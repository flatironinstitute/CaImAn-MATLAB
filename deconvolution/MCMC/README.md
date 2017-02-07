[![Join the chat at https://gitter.im/epnev/ca_source_extraction](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/epnev/ca_source_extraction?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

MCMC spike inference in continuous time
==========================

The code takes as an input a time series vector of calcium observations
and produces samples from the posterior distribution of the underlying
spike in continuous time. The code also samples the model parameters
(baseline, spike amplitude, initial calcium concentration, firing rate,
noise variance) and also iteratively re-estimates the discrete time
constant of the model. More info can be found at

Pnevmatikakis, E., Merel, J., Pakman, A. &amp; Paninski, L. (2014).
Bayesian spike inference from calcium imaging data. Asilomar Conf. on
Signals, Systems, and Computers. http://arxiv.org/abs/1311.6864

For initializing the MCMC sampler, the algorithm uses the constrained deconvolution method maintained separately in https://github.com/epnev/constrained-foopsi

### Contributors

Eftychios A. Pnevmatikakis, Simons Foundation

John Merel, Columbia University

### Contact

For questions join the chat room (see the button above) or open an issue (for bugs etc).

### Acknowledgements

Special thanks to Tim Machado for providing the demo dataset.

License
=======

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
