uniform mat4 modelview;
attribute vec4 pos;
attribute vec2 vtex;
varying vec2 ftex;

void main()
{

    gl_Position = modelview * pos;
    ftex = vtex;
}
