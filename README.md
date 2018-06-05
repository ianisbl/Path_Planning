# Path_Planning
by Ianis Bougdal-Lamber and Gael Colas graduate students at Stanford.

This repository gathers our final projects for the AA222: "Engineering Design Optimization" and the AA203 : "Optimal Control" classes at Stanford School of Engineering (2018). Our professors were Pr. Mykel Kochenderfer and Pr. Marco Pavone.

Languages: Python, Matlab

Goal: solve the 2D optimal Path Planning problem using different optimization methods.
This problem can be defined as follows: find the shortest feasible path between an initial and a final cell on a 2D map.
A feasible path does not collide with obstacles. 

This project involved 2 distinct parts:
  - Genetic Algorithm implementation and analysis: the AA222 project used a Genetic Algorithm to evolve a population of path individuals into an optimal path. A multi-objectve version was also considered where the goal is to minimize both the length and the difficulty of the path.
  - Differential Dynamic Programming implementation and analysis:

Our algorithms were compared to the state-of-the-art: the famous A* algorithm.
