## Animating a Tippe Top using Numerical Methods

This project is a study of the dynamics of a special kind of spinning top known as a tippe top. A tippe top has a particular mass distribution such that it initially spins on its spherical end, then while spinning, the tippe top inverts itself and continues spinning on its extended stem. See this YouTube video for an example: https://youtu.be/AyAgeUneFds

My full paper is available: Joneja-TippeTop.docx or .pdf
The slides from my final presentation on May 8, 2017: Joneja-TippeTop-final.pptx

### Simple spinning top
A model of a simple, symmetric spinning top was explored first. The code structure involves solving a system of ODEs using fixed step size RK4, transforming the solution from the Euler angle coordinate system to the 'lab' reference frame. Animations are produced using VPython.

Click any of the following to watch the animation results on YouTube.

[![youtube](http://img.youtube.com/vi/I99i_cIsqz0/0.jpg)](https://www.youtube.com/watch?v=I99i_cIsqz0)
[![youtube](http://img.youtube.com/vi/35D7rX08HPk/0.jpg)](https://www.youtube.com/watch?v=35D7rX08HPk)
[![youtube](http://img.youtube.com/vi/u3LqBxvr-7o/0.jpg)](https://www.youtube.com/watch?v=u3LqBxvr-7o)

simpletop.py produces plots of nutation, rate of precession and spin rate over time. In addition, plots of Energy over time are also produced.
simpletop_v.py runs the animation (provided you have Python 2.7 with the visual package installed).

### Tippe top
Next, the full model of the tippe was studied. The code structure is similar to that of the simple spinning top presented above. 
The key differences are the model of the tippe top itself and the fact that an adaptive RK4 method was used here instead. 

Click any of the following to watch the animation results on YouTube.

[![youtube](http://img.youtube.com/vi/ndNa3KMrU84/0.jpg)](https://www.youtube.com/watch?v=ndNa3KMrU84)
[![youtube](http://img.youtube.com/vi/fjGDBIcBP9U/0.jpg)](https://www.youtube.com/watch?v=fjGDBIcBP9U)
[![youtube](http://img.youtube.com/vi/7pJ0hYqawHk/0.jpg)](https://www.youtube.com/watch?v=7pJ0hYqawHk)

tippe_adaptive.py produces plots of nutation, rate of precession and spin rate over time. 
tippe_adaptive_v.py runs the animation (again, Python 2.7 with visual package required).