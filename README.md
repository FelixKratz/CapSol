# CapSol
## Compile the CapSol Python module
### macOS
* Required dependencies ([brew](https://brew.sh) needs to be installed)
```bash
brew install gsl
brew install armadillo
brew install pybind11
```
* Compile the code with `make`. This will put a compiled static object file (`.so`) in the `frontend` folder. Every python script using the CapSol module needs to have this object file in its search path (e.g. by residing in the same folder).

### Linux
* Install `gsl`, `armadillo` and `pybind11` with the package manager of choice.
* Compile the code with `make`. This will put a compiled static object file (`.so`) in the `frontend` folder. Every python script using the CapSol module needs to have this object file in its search path (e.g. by residing in the same folder).

### Testing the Python module
If you want to run a functionality test of the module you can simply run `make clean && make test`. This will run a short check of the most relevant systems and report them in the commandline.

## Using the CapSol Python module
Some examples are provided in the `frontend` folder of this repository as Jupyter Notebooks (interactive python runtime).

### Importing the module
The module can be imported into a python script (which has the static object file `.so` compiled above in its search path) simply by calling:
```python
from CapSol import <module>
```
where `<module>` is any of the following available submodules of the package.

## Available `<module>`s
### The `Laplace` module
The `Laplace` module makes available the numerical solver for the Young-Laplace
equation together with the numerical fitting routine for these shape equations:

* Generating shapes for given parameters `p_L`, `rho` and `omega` (these are exactly $\tilde{p}_L, \Delta\tilde{\rho}$ and $\Omega$ from Chap. 1 of the Dissertation)
```python
laplace = Laplace(p_L, rho, omega)
```

* Fitting parameters for a given `shape`, consisting of an array of non-dimensional x-y shape coordinate pairs (the `shape` array can be generated from images by using the provided `image` package detailed above)
```python
laplace = Laplace(shape)
```

The properties of the resulting `laplace` object are:
| Property | Description |
| -------- | ----------- |
| `laplace.p_L` | The dimensionless apex pressure |
| `laplace.rho` | The dimensionless density contrast |
| `laplace.shape` | An array of shape coordinates |
| `laplace.valid` | Boolean indicating if the solution is valid |

### The `Hooke` module
The `Hooke` module makes available the numerical solver for the non-linear Hookean shape
equations together with the numerical fitting routine for these shape equations:

* Generating shapes for given parameters `p_L`, `rho`, `K`, `nu`, `tau_s_0` (these are exactly $\tilde{p}_L$, $\Delta\tilde{\rho}$, $K_{2D} / \gamma$, $\nu_{2D}$ and $\tau_s(s = 0) / \gamma$ from Chap. 1 of the Dissertation)
```python
hooke = Hooke(K, nu, tau_s_0, p_L, rho)
```

* Fitting parameters for a given `shape` and given reference shape parameters `p_L` and `rho` (to be obtained e.g. by the `Laplace` module)
```python
hooke = Hooke(shape, p_L, rho)
```

The properties of the resulting `hooke` object are:
| Property | Description |
| -------- | ----------- |
| `hooke.p_L` | The dimensionless apex pressure of the reference shape |
| `hooke.p_a` | The dimensionless apex pressure of the deformed shape |
| `hooke.rho` | The dimensionless density contrast |
| `hooke.K` | The dimensionless two dimensional compression modulus |
| `hooke.nu` | The dimensionless two dimensional Poisson's ratio |
| `hooke.tau_s_0` | The dimensionless apex stress |
| `hooke.tau_s` | An array of dimensionless meridional stresses along the shape |
| `hooke.tau_phi` | An array of dimensionless circumferential stresses along the shape |
| `hooke.lambda_s` | An array of meridional stretches along the shape |
| `hooke.lambda_phi` | An array of circumferential stretches along the shape |
| `hooke.shape` | An array of shape coordinates of the deformed shape |
| `hooke.valid` | Boolean indicating if the solution is valid |

### The `Kelvin` module
The `Kelvin` module makes available the numerical solver for the viscoelastic Kelvin-Voigt shape
equations together with the numerical fitting routine for these shape equations:

* Generating shapes for given parameters `p_L`, `rho`, `K`, `nu`, `eta`, `apexStresses`, `periods` (these are exactly $\tilde{p}_L$, $\Delta\tilde{\rho}$, $K_{2D} / \gamma$, $\nu_{2D}$ and $\eta_{2D} \omega / \gamma$ together with an array of dimensionless apex stresses and the number of oscillation periods (see Chap. 1 of the Dissertation))
```python
kelvin = Kelvin(apexStresses, K, nu, eta, p_L, rho, periods)
```

* Fitting parameters for a given `shapes` (array of `shape` used above), given reference shape parameters `p_L` and `rho` (to be obtained e.g. by the `Laplace` module), and given number of oscillation `periods`
```python
kelvin = Kelvin(shapes, p_L, rho)
```

The properties of the resulting `kelvin` object are:
| Property | Description |
| -------- | ----------- |
| `kelvin.p_L` | The dimensionless apex pressure of the reference shape |
| `kelvin.rho` | The dimensionless density contrast |
| `kelvin.K` | The dimensionless two dimensional compression modulus |
| `kelvin.nu` | The dimensionless two dimensional Poisson's ratio |
| `kelvin.eta` | The dimensionless viscosity |
| `kelvin.apexStresses` | An array of the non-dimensional apex stresses of the deformed shapes |
| `kelvin.shapes` | An array of shape coordinates of the deformed shape |
| `kelvin.valid` | Boolean indicating if the solution is valid |

### The `Contact` module
The `Contact` module makes available the numerical solver for the contact problem discussed in Chap. 4 of the Dissertation

* Generating shapes for given parameters `p_L_u`, `rho_u`, `p_L_d`, `rho_d`, `K_u`, `nu_u`, `tau_s_0_u`, `K_d`, `nu_d`, `tau_s_0_d`, `contact_length`, `tension_ratio`, `tension_inner_u`, `tension_inner_d`, `tension_ud`  (these are exactly $\tilde{p}_L^u$, $\Delta\tilde{\rho}^u$, $\tilde{p}_L^d$, $\Delta\tilde{\rho}^d$, $K_{2D}^u / \gamma^u$ , $\nu_{2D}^u$, $\tau_s^u(s = 0) / \gamma^u$, $K_{2D}^d / \gamma^u$ , $\nu_{2D}^d$, $\tau_s^d(s = 0) / \gamma^u$, l, $\Gamma$, $\gamma^u_i / \gamma^u$,  $\gamma^d_i / \gamma^u$ and $\gamma^{ud} / \gamma^u$, see Chap. 4 of the Dissertation)
```python
contact = Contact(p_L_u, rho_u, p_L_d, rho_d,
                  nu_u, K_u, tau_s_0_u,
                  nu_d, K_d, tau_s_0_d,
                  contact_length, tension_ratio,
                  tension_inner_u, tension_inner_d,
                  tension_ud)
```

The properties of the resulting `contact` object are:
| Property | Description |
| -------- | ----------- |
| `contact.p_a_u` | The dimensionless apex pressure of the upper deformed shape |
| `contact.p_a_d` | The dimensionless apex pressure of the lower deformed shape |
| `contact.tau_s_u` | An array of dimensionless meridional stresses along the upper shape |
| `contact.tau_s_d` | An array of dimensionless meridional stresses along the lower shape |
| `contact.tau_s_u` | An array of dimensionless circumferential stresses along the upper shape |
| `contact.tau_s_d` | An array of dimensionless circumferential stresses along the lower shape |
| `contact.lambda_s_u` | An array of meridional stretches along the upper shape |
| `contact.lambda_s_d` | An array of meridional stretches along the lower shape |
| `contact.lambda_phi_u` | An array of circumferential stretches along the upper shape |
| `contact.lambda_phi_d` | An array of circumferential stretches along the lower shape |
| `contact.r_u` | An array of radial shape coordinates of the upper deformed shape |
| `contact.r_d` | An array of radial shape coordinates of the lower deformed shape |
| `contact.z_u` | An array of height shape coordinates of the upper deformed shape |
| `contact.z_d` | An array of height shape coordinates of the lower deformed shape |
| `contact.f` | The dimensionless contact force |
| `contact.valid` | Boolean indicating if the solution is valid |

## Helper module for image processing
If you want to process images into a `shape` representation you can use the provided image processing routine contained in the `frontend` folder like so:
```python
from image import loadImage, loadImages, Mask
```
Then you will be able to load an image with a given file `path` like this:
```python
mask = Mask(100, 100, 0, 100)
shape = loadImage(path, flip=True, mask=mask, verbose=True)
```
where you can mask the image with the provided `Mask` class to exclude certain regions of the image. You can flip the image with the `flip` property (to be used when experimenting in a rising geometry) and show the loaded, masked and preprocessed image with `verbose=True` (see for example `frontend/LaplaceFit.ipynb`).

To load an array of images in parallel from an array of file `paths` you can use:
```python
shapes = loadImages(paths, flip=True, mask=mask)
```

