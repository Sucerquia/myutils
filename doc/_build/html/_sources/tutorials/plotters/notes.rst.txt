Standard
========

Data
----

one list (array)
****************

If one list (or array) is the only argument, it is assumed that corresponds with
the data in the y axis corresponding with the number label in the x axis. Then,
in the next example, we are going to generate 200 random numbers and to pass
that list of data as the only argument. As you can see, the y axis corresponds
with the generated random numbers and the x axis goes from 1 to 200.

.. literalinclude:: one_array.py

Note: You can also concatenate several lists in one list. But we highly
recommend to separate each of then in case it is practical. It allows you to
have more control in your plot.

two lists (arrays)
******************

In case of two list, it is assumed that first of them corresponds with the list
of x-coordinates and the other corresponds with the list of y-coordinates. Of
course, those lists must have the same lenght.

As an example, let's take a set of 200 random numbers and create a growing list
of numbers by doing a cummulative sum. Then, let's fit this data into an
hypertangent function composed with a linear function with m=1/25 and b=-2.

.. literalinclude:: hypertangent.py

Adding several curves to the same plot
**************************************

A big advantage of plots is to compare the behavior of two set of data and the
best way to do it is to put them at the same plot. To do so using standard, you
only have to call 
2COMPLETE add link to function standard
standard several times. 