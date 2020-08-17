# Computational Physics

The material here aims at giving you an introduction to several of the most used algorithms in Computational Science. These algorithms cover topics such as advanced numerical integration using Gaussian quadrature, Monte Carlo methods with applications to random processes, Markov chains, integration of multidimensional integrals and applications to problems in statistical physics and quantum mechanics. Other methods which are presented are eigenvalue problems, from the simple Jacobi method to iterative Krylov methods. Popular methods from linear algebra are also discussed. A good fraction of the course is also devoted to solving ordinary differential equations with or without boundary conditions and finally methods for solving partial differential equations. You will also find material on popular Machine Learning algorithms, starting with various linear regression methods and ending with neural networks. The focus for the Machine Learning algorithms is on supervised learning.

The course is project based and through various projects, normally four to five, you will be exposed to fundamental research problems from various fields (Physics, Geophysics, Chemistry, Mathematics, Statistics etc), where, if possible, we aim at reproducing state of the art scientific results. You will learn to develop and structure codes when solving the projects, develop a critical understanding of the strengths and limits of the various numerical methods, become familiar with supercomputing facilities and parallel computing and learn to write scientific projects. 

## Instructors information
* _Name_: Morten Hjorth-Jensen
  * _Email_: morten.hjorth-jensen@fys.uio.no
  * _Phone_: +47-48257387
  * _Office_: Department of Physics, University of Oslo, Eastern wing room 470 
  * _Office hours_: *Anytime*! In Fall Semester 2020 (FS20), as a rule of thumb office hours are planned via computer or telephone. Individual or group office   hours will be performed via zoom. Feel free to send an email for planning. In person meetings may also be possible if allowed by the University of Oslo's COVID-19 instructions (see below for links).
* _Name_: Anders Kvellestad
  * _Email_: anders.kvellestad@fys.uio.no
  * _Office_: Department of Physics, University of Oslo, Eastern wing room 447 

##  Teaching Assistants FS20
* Sebastian Wither-Larsen, sebastwi@student.matnat.uio.no
  * _Office_: Department of Physics, University of Oslo, Eastern wing room 454
* René Alexander Ask, r.a.ask@fys.uio.no
* Kaspara Gåsvær, k.s.gasvar@fys.uio.no
* Maria Linea Horgen, m.l.horgen@fys.uio.no
* Aksel Graneng, akselgraneng@gmail.com

## Practicalities
This course will be delivered in a hybrid mode, with online lectures and on site or online laboratory sessions. 

1. Four lectures per week, Fall semester, 10 ECTS. The lectures will be fully online. The lectures will be recorded and linked to this site and the official University of Oslo website for the course;
2. Two hours of laboratory sessions for work on computational projects for each group. Due to social distancing, at most 15 participants can attend. There will  also be fully digital laboratory sessions for those who cannot attend physically;
3. Three projects which are graded and count 1/3 each of the final grade, five projects in total;
4. The course is offered as  FYS4150 (Master of Science level) and  FYS3150 (senior undergraduate level);
5. We use Piazza for course communication, a special link on how to register to Piazza can be found at the official University of Oslo page for the course or just use the link here https://piazza.com/uio.no/fall2020/fys3150. Slack is also used for course communication. The Slack link is https://compphysicsuio.slack.com ;
6. Videos of teaching material are available via the links at https://compphysics.github.io/MachineLearning/doc/web/course.html;
7. Weekly emails with summary of activities will be mailed to all participants;

## Grading
Grading scale: Grades are awarded on a scale from A to F, where A is the best grade and F is a fail. There are three projects which are graded and each project counts 1/3 of the final grade. The total score is thus the average from all three projects.

The final number of points is based on the average of all projects (including eventual additional points) and the grade follows the following table:

 * 92-100 points: A
 * 77-91 points: B
 * 58-76 points: C
 * 46-57 points: D
 * 40-45 points: E
 * 0-39 points: F-failed

## Required Technologies

Course participants are expected to have their own laptops/PCs. We use _Git_ as version control software and the usage of providers like _GitHub_, _GitLab_ or similar are strongly recommended.

We will make extensive use of C++ and/or Python as programming language and its
myriad of available libraries.  You will find
Jupyter notebooks invaluable in your work.  You can run _R_
codes in the Jupyter/IPython notebooks, with the immediate benefit of
visualizing your data. You can also use compiled languages like C++,
Rust, Julia, Fortran etc if you prefer. 


If you have Python installed (we strongly recommend Python3) and you feel
pretty familiar with installing different packages, we recommend that
you install the following Python packages via _pip_ as 

* pip install numpy scipy matplotlib ipython scikit-learn mglearn sympy pandas pillow 

For OSX users we recommend, after having installed Xcode, to
install _brew_. Brew allows for a seamless installation of additional
software via for example 

* brew install python3

For Linux users, with its variety of distributions like for example the widely popular Ubuntu distribution,
you can use _pip_ as well and simply install Python as 

* sudo apt-get install python3

### Python installers

If you don't want to perform these operations separately and venture
into the hassle of exploring how to set up dependencies and paths, we
recommend two widely used distrubutions which set up all relevant
dependencies for Python, namely 

* Anaconda:https://docs.anaconda.com/, 

which is an open source
distribution of the Python and R programming languages for large-scale
data processing, predictive analytics, and scientific computing, that
aims to simplify package management and deployment. Package versions
are managed by the package management system _conda_. 

* Enthought canopy:https://www.enthought.com/product/canopy/ 

is a Python
distribution for scientific and analytic computing distribution and
analysis environment, available for free and under a commercial
license.

Furthermore, Google's Colab:https://colab.research.google.com/notebooks/welcome.ipynb is a free Jupyter notebook environment that requires 
no setup and runs entirely in the cloud. Try it out!

### Useful Python libraries
Here we list several useful Python libraries we strongly recommend (if you use anaconda many of these are already there)

* _NumPy_:https://www.numpy.org/ is a highly popular library for large, multi-dimensional arrays and matrices, along with a large collection of high-level mathematical functions to operate on these arrays
* _The pandas_:https://pandas.pydata.org/ library provides high-performance, easy-to-use data structures and data analysis tools 
* _Xarray_:http://xarray.pydata.org/en/stable/ is a Python package that makes working with labelled multi-dimensional arrays simple, efficient, and fun!
* _Scipy_:https://www.scipy.org/ (pronounced “Sigh Pie”) is a Python-based ecosystem of open-source software for mathematics, science, and engineering. 
* _Matplotlib_:https://matplotlib.org/ is a Python 2D plotting library which produces publication quality figures in a variety of hardcopy formats and interactive environments across platforms.
* _Autograd_:https://github.com/HIPS/autograd can automatically differentiate native Python and Numpy code. It can handle a large subset of Python's features, including loops, ifs, recursion and closures, and it can even take derivatives of derivatives of derivatives
* _SymPy_:https://www.sympy.org/en/index.html is a Python library for symbolic mathematics. 
* _scikit-learn_:https://scikit-learn.org/stable/ has simple and efficient tools for machine learning, data mining and data analysis
* _TensorFlow_:https://www.tensorflow.org/ is a Python library for fast numerical computing created and released by Google
* _Keras_:https://keras.io/ is a high-level neural networks API, written in Python and capable of running on top of TensorFlow, CNTK, or Theano
* And many more such as _pytorch_:https://pytorch.org/,  _Theano_:https://pypi.org/project/Theano/ etc 


## Personal Hygiene
All participants attending the laboratory sessions must maintain proper hygiene and health practices, including:
* frequently wash with soap and water or, if soap is unavailable, using hand sanitizer with at least 60% alcohol;
* Routinely cleaning and sanitizing living spaces and/or workspace;
* Using the bend of the elbow or shoulder to shield a cough or sneeze;
* Refraining from shaking hands;

## Adherence to Signage and Instructions 
Course participants  will (a) look for instructional signs posted by UiO or public health authorities, (b) observe instructions from UiO or public health authorities that are emailed to my “uio.no” account, and (c) follow those instructions.
The relevant links are https://www.uio.no/om/hms/korona/index.html and https://www.uio.no/om/hms/korona/retningslinjer/veileder-smittevern.html

## Self-Monitoring
Students will self-monitor for flu-like symptoms (for example, cough, shortness of breath, difficulty breathing, fever, sore throat or loss of taste or smell). If a student experiences any flu-like symptoms, they will stay home and contact a health care provider to determine what steps should be taken.
## Exposure to COVID-19 
If a student is exposed to someone who is ill or has tested positive for the COVID-19 virus, they will stay home, contact a health care provider and follow all public health recommendations. You may also contact the study administration of the department where you are registered as student. 
## Compliance and reporting 
Those who come to UiO facilities must commit to the personal responsibility necessary for us to remain as safe as possible, including following the specific guidelines outlined in this syllabus and provided by UiO more broadly (see links below here). 

## Additional information
See https://www.uio.no/om/hms/korona/index.html and https://www.uio.no/om/hms/korona/retningslinjer/veileder-smittevern.html. For English version, click on the relevant link.
