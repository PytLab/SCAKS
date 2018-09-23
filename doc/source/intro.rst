.. currentmodule:: MiKiAC

Introduction
=============

Welcome to the `MiKiAC` documentation!

`MiKiAC` is a
Python implementation of the Micro-kinetic model solving model with a user-friendly 
web GUI distributed under the GPLv3 license. It's an acronym for Micro Kinetics 
Analysis for Catalyst.

You can always find the latest stable version of the program
here: https://github.com/pytlab/mikiac

This documentation describes version 1.0.0


Why Micro-kinetic Model for Catalyst?
-------------------------------------
Rational design of catalysis assisted by means of first-principle calculation is
current one of the most important topic in field of heterogeneous catalysis 
because the traditional try-and-error method can’t meet the rapidly increasing 
demands of catalyst industry development. In the process of rational design and 
catalysis screening implementation, solving the micro-kinetic model provides the 
theoretical basis for describing the turnover frequency and selectivity.

Why `MiKiAC`?
---------------------------
Despite the great development of microkinetics simulation methods and software, 
the implementations of kinetics analysis with relatively fixed code architect are 
not flexible and scalable enough for more and more complex catalytic system simulation. 
Besides, the absence of user-friendly interfaces is also the obstacle on the road 
for chemical researchers without programming knowledges to use those programs. 

To this end, we present a Python module called “MiKiAC” with web GUI to help 
researchers to solve microkinetics more easily by lowering the learning and using 
barriers. As a complete Python module which can be imported in other user customized 
programs, MiKiAC provide robust and flexible interfaces to help expert users 
create one or more models and solve them at a time with many components built in 
such as powerful reaction expression parser, object-oriented based energy profile 
plotter and so on. In order to provide a more friendly and interactive user 
interfaces for other chemical researchers, we also use the famous Python micro 
web framework Flask and the web front-end framework Bootstrap to build a web 
application using MiKiAC as the back-end calculation core. Then we can run the 
microkinetic model application on both local and remote server and make it possible 
for a scalable and large scaled cloud computing service.

Other implementation codes
...........................
If you don't find `MiKiAC` suitable for your purpose, I strongly ecourage you to 
consider any of the other available micro-kinetic model solving codes such as 
`CatMAP`_, `CHEMKIN`_ and so on before implementing your own methods.

.. _CatMAP: https://github.com/SUNCAT-Center/catmap
.. _CHEMKIN: http://www.reactiondesign.com/products/chemkin

