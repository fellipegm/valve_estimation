# Valve Estimation

This repository is the porting of my PhD Thesis from MATLAB to C++ to obtain higher computational performance.

## Instalation

I used the Visual Studio Community 2019 with Intel Compiler 19.1 for Windows or just the Intel Compiler 19.1 in Linux. After some changes in the code, I think that it only works in Windows.

It seems that it will not take too much effort to make it work in Linux again, mainly I made checks relying in windows headers to verify if the passed path exists.

### Usage

To change the parameters of the simulations or estimations, it is necessary to recompile the program. This may be updated in future versions.

The program is designed to work with CLI.

A help is displayed with ValveEstimation.exe or ValveEstimation.exe --help

```bash
> ValveEstimation.exe # or ValveEstimation.exe --help

Usage: Usage: ValveEstimation --directory=<path> --type=<type> --excitation=<excitation> --models=<model1-model2-etc> -k --noise -n 1 --valve=<valve> --load=<real observations csv file with path>
User must define the type of experiment with --type=<type>
<type> can be: simulation       estimation      real-estimation real-simulation cl-simulation   cl-estimation

User must define the type of excitation signal with --excitation=<excitation>
<excitation> can be: sinusoidal aleatory

User must define the friction models separated with minus with --models=<model1-model2>
models can be: kano     he      choudhury       karnopp lugre   gms     sgms    gms1

-k is used to estimate k and F_init toghether
--noise is used to simulate noise in the estimations
-n is used to set how many estimations are necessary -n 10
--valve=<graphite,teflon> sets the valve data used in real estimation or simulation
--load=<file with path> sets the file with real data for estimation
```

If the user wants to simulate the he and kano models (using the compiled parameters), saving the output files in D:\tests folder and using the sinusoidal excitation signal, the command is:

```bash
> ValveEstimation.exe --directory=D:\tests\ --type=simulation --excitation=sinusoidal --models=he-kano
```

## Work in Progress

* [ ] Closed loop algorithms
  * [ ] Algorithms to simulate the valve in closed loop
  * [ ] Algorithms to obtain the initial time where the excitation started
* [ ] Load the valve data from an external file

* [ ] Describe the expected data format

* [ ] Use binary data to avoid spending too much space in hdd

## License

MIT
