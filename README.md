# SpaceGame
A Space Game made with Unity.

<img src="https://github.com/bogger12/SpaceGame/blob/3b6e416339a3325e6d3ccc5b1592ef55c65ba322/example.png" width="500" height="500*1.14"></img>
## To Do List
- [ ] Add (simple) atmosphere sprite to planet
- [ ] Add apsis markers
- [x] Add zoom in and out.
- [x] Add moon
- [x] Add ability to focus on different celestial objects
- [x] Add hyperbolic and parabolic trajectories
- [x] Add SOI of celestial objects such as moons, and get trajectory past them
- [ ] Add manuever planners
- [x] Fix crash when orbit goes to 0 width
- [ ] Add proper controls - z,x, others
- [ ] Add UI for stuff like fuel, crew, thrust, position, orientation, stats
- [ ] Make planets have dynamic lighting depending on where they are facing (shaders?)
- [x] Add sun and more planets in system

### Future Ideas
- [ ] Add higher res pngs (svgs maybe) for planets, so we can land on them etc
- [ ] Improve rocket sprite
- [ ] Add animations for rocket etc. (particles?)
- [ ] Add crew and EVA capabilities

### Potential Ideas
- [ ] Add refueling capabilities on planets
- [ ] Add more to do in eva, such as gathering materials, npcs
- [ ] Add space wreckage
- [ ] Add some story maybe? Motivation to go explore


## Resource List

**[Astrophysics Equations](https://static1.squarespace.com/static/54b38552e4b055a31e5e3e47/t/5b5c071970a6addd33a8a85c/1532757806053/ultimate-astrophysics-cheat-sheet_1-0.pdf)**\
[Orbital Mechanics Reference](https://orbital-mechanics.space/intro.html)
___
### General

<details open>
  
**<summary>Orbits/Gravity</summary>**
  
* https://en.wikipedia.org/wiki/Earth
* https://en.wikipedia.org/wiki/Gravitational_constant
* https://en.wikipedia.org/wiki/Orbit_equation
* https://en.wikipedia.org/wiki/Hyperbolic_trajectory
* https://en.wikipedia.org/wiki/Parabolic_trajectory
* https://en.wikipedia.org/wiki/Two-body_problem
* https://en.wikipedia.org/wiki/Elliptic_orbit
* https://en.wikipedia.org/wiki/Ellipse
* https://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion
* https://en.wikipedia.org/wiki/Orbital_elements
* https://en.wikipedia.org/wiki/Oberth_effect
* https://en.wikipedia.org/wiki/Sphere_of_influence_(astrodynamics)
* https://en.wikipedia.org/wiki/Hill_sphere
* https://en.wikipedia.org/wiki/Roche_limit
* https://en.wikipedia.org/wiki/Titius%E2%80%93Bode_law
</details>

<details open>
  
**<summary>Vectors</summary>**

* [Vector Definition](https://en.wikipedia.org/wiki/Vector_(mathematics_and_physics))
* [Unit Vector](https://en.wikipedia.org/wiki/Unit_vector)
* [Cross Product](https://en.wikipedia.org/wiki/Cross_product)
* [Dot Product](https://en.wikipedia.org/wiki/Dot_product)
* **[Vector Lesson](http://physics.bu.edu/~redner/211-sp06/class03/comp_vectors.html)**
</details>


___
### Orbital Elements:
* $e$&nbsp;&nbsp;[Eccentricity](https://en.wikipedia.org/wiki/Orbital_eccentricity) | $\textbf{e}$&nbsp;&nbsp;[Eccentricity Vector](https://en.wikipedia.org/wiki/Eccentricity_vector)
* $a$&nbsp;&nbsp;[Semi Major Axis](https://en.wikipedia.org/wiki/Semi-major_and_semi-minor_axes)
* $i$&nbsp;&nbsp;[Inclination](https://en.wikipedia.org/wiki/Orbital_inclination)
* $Ω$&nbsp;&nbsp;[Longtitude of the Ascending Node](https://en.wikipedia.org/wiki/Longitude_of_the_ascending_node)
* $ω$&nbsp;&nbsp;[Argument of Periapsis](https://en.wikipedia.org/wiki/Argument_of_periapsis)
* $ν$&nbsp;&nbsp;[True Anomaly](https://en.wikipedia.org/wiki/True_anomaly) (usually at epoch)

###### Other Important Elements:
* $t_0$&nbsp;&nbsp;[Epoch](https://en.wikipedia.org/wiki/Epoch_(astronomy))
* $T$&nbsp;&nbsp;[Orbital Period](https://en.wikipedia.org/wiki/Orbital_period)
* $M$&nbsp;&nbsp;[Mean Anomaly](https://en.wikipedia.org/wiki/Mean_anomaly)
* $E$&nbsp;&nbsp;[Eccentric Anomaly](https://en.wikipedia.org/wiki/Eccentric_anomaly)
* $h$&nbsp;&nbsp;[Specific Angular momentum](https://en.wikipedia.org/wiki/Specific_angular_momentum)

___

### Helpful Articles about Programming:
* https://mikhail-szugalew.medium.com/simulating-gravity-in-unity-ae8258a80b6d
* https://phys.libretexts.org/Bookshelves/Astronomy__Cosmology/Celestial_Mechanics_(Tatum)/09%3A_The_Two_Body_Problem_in_Two_Dimensions/9.08%3A_Orbital_Elements_and_Velocity_Vector#mjx-eqn-9.5.31
* https://space.stackexchange.com/questions/1904/how-to-programmatically-calculate-orbital-elements-using-position-velocity-vecto
* https://space.stackexchange.com/questions/52090/how-can-i-calculate-the-future-position-of-a-satellite-orbiting-a-central-body-a
* https://space.stackexchange.com/questions/4093/how-can-i-predict-what-an-objects-orbital-state-vectors-will-be-in-the-future?rq=1
* https://space.stackexchange.com/questions/3781/calculating-the-time-since-periapse-of-an-orbital-position
* https://physics.stackexchange.com/questions/333876/how-does-one-calculate-the-time-to-go-from-one-point-in-an-orbit-to-another#:~:text=tp(r)%3D√,tp(t1)
* **https://space.stackexchange.com/questions/8911/determining-orbital-position-at-a-future-point-in-time**

Orbital Elements Calculator:
* http://orbitsimulator.com/formulas/OrbitalElements.html

This Lecture I found which has lots of useful stuff:
* https://control.asu.edu/Classes/MAE462/462Lecture05.pdf

Orbital Mechanics Educational Website
* https://orbital-mechanics.space/intro.html