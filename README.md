# occam1dmt
Toy Occam inversion code for 1D MT (Magnetotellurics) method in Matlab

## DATA
The example data are (literally) copied from Table 5 of the Constable 1987 paper, which is in turn from Cull (1985). See: 

Constable, S. C., Parker, R. L., & Constable, C. G. (1987). Occam’s inversion: A practical algorithm for generating smooth models from electromagnetic sounding data. Geophysics, 52(3), 289–300. 

Cull, J. P. (1985). Magnetotelluric soundings over a Precambrian contact in Australia. Geophys. J. Roy. Astr. Sot., 80, 661-675.

see also my toy occam code for DC resistivity:

[https://github.com/dong-hao/occam1ddc]

## USAGE
See example/testbench.m for a simple demonstration on how to load the data and call the inversion code. 

## something like a disclaimer

This was one of many toy codes I fiddled with when I was a student - I hope this could be useful to our students nowadays in the EM community. 
Those who want to try this script are free to use it on academic/educational cases. But of course, I cannot guarantee the script to be working properly and calculating correctly (although I wish so). Have you any questions or suggestions, please feel free to contact me (but don't you expect that I will reply quickly!).  

## HOW TO GET IT
```
git clone https://github.com/dong-hao/occam1dmt/ your_local_folder
```

## UNITS
The internal scale here is log10(Ohmm) - instead of linear scale for both resistivity and apparent resistivity. The layer depth is in (linear) metres, while the impedance phase is in rads. 

## ERRORS    
Currently the internal error here is standard deviation.

## HOW TO GET UPDATED
```
cd to_you_local_folder
git pull 
```

## Contact

DONG Hao –  donghao@cugb.edu.cn

China University of Geosciences, Beijing 

Distributed under the GPL v3 license. See ``LICENSE`` for more information.

[https://github.com/dong-hao/occam1dmt]

## Contributing

Those who are willing to contribute are welcomed to try - but I probably won't have the time to review the commits frequently (not that I would expect there will be any). 

1. Fork it (<https://github.com/dong-hao/occam1dmt/fork>)
2. Create your feature branch (`git checkout -b feature/somename`)
3. Commit your changes (`git commit -am 'Add some features'`)
4. Push to the branch (`git push origin feature/somename`)
5. Create a new Pull Request - lather, rinse, repeat 
