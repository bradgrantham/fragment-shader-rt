Here's a stab at usable source code for my fragment-program sphere
ray-tracer, originally posted at http://plunk.org/~grantham/ray1 .

The only currently working targets for the Makefile are "make_js_world"
and "web_models".

To make a working version of the ray-tracer web page, do the following:
* clone this repository
* make web_models
* python -m SimpleHTTPServer
* point your WebGL-capable browser at the address printed
