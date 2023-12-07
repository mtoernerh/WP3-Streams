# WP3-Streams (WaterPlan 3)
This is a an old project for including +3000 km streams targeted by the Water Framework Directive (WFD) to the Danish National Water Resources Model (DK-model). The state of the code is chaotic, but a working example for the island of Fyn is included. The objective of the code is to process shapefile linestrings representing streams, but that are lacking proper stream delineation and flow direction, and correcting those two attributes. Furthermore, the function attempts to resolve looping- and self-intersecting streams although proper delination and flow direction is ambigious in those cases.

It is important that the linestring input is topologically connected, that is the function will not be able to detect upstream segments unless they are actually connected. This could perhaps be somewhat circumvented by using a buffer search around the start vertex, and snapping to the potential nearest end vertices of the upstream linestrings.


![Description](images/WP3_description.svg)

Below is a simple example of processing streams represented as simple linestrings, where proper delineation of main- and substreams and flow direction is assigned. Each color represents a unique linestring, and the mid-line arrow shows the direction of flow. You will notice in the before image, that flow direction is random with no apparent pattern, furthermore the lack of stream delineation makes it ambigious as to which linestrings belong to which streams. This if of course one intepretation of that, but it is easily understandable and the method for delineating is systematic and consistent.



![Example](images/stream_correction.svg)








[Example]: https://github.com/mtoernerh/WP3-Streams/images/stream_correction.svg
