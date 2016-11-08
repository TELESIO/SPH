#include<iostream>
#include <stdio.h>
#include <stdlib.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <time.h>
#include "engine.h"
#include<memory>

#define DT (0.001d)
#define WX (0.5d)
#define WY (0.5d)
#define WZ (0.5d)

const double MASS=0.001;
const double PRADIUS=0.0157;


struct ModelView{
    GLfloat x_rot;
    GLfloat y_rot;
    GLfloat z_trans;
};

struct ModelView model_view;
time_t start_time, end_time;

SPHEngine* engine;


/* Axes display list */
static GLuint axes_list;

void initAxesList(){

    /* Create a display list for drawing axes */
    axes_list = glGenLists(1);
    glColor3b(1,1,1);
    glNewList(axes_list, GL_COMPILE);

    glColor4ub(0, 0, 255, 255);
    glBegin(GL_LINE_STRIP);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(1.0f, 0.0f, 0.0f);
    glVertex3f(0.75f, 0.25f, 0.0f);
    glVertex3f(0.75f, -0.25f, 0.0f);
    glVertex3f(1.0f, 0.0f, 0.0f);
    glVertex3f(0.75f, 0.0f, 0.25f);
    glVertex3f(0.75f, 0.0f, -0.25f);
    glVertex3f(1.0f, 0.0f, 0.0f);
    glEnd();
    glBegin(GL_LINE_STRIP);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 1.0f, 0.0f);
    glVertex3f(0.0f, 0.75f, 0.25f);
    glVertex3f(0.0f, 0.75f, -0.25f);
    glVertex3f(0.0f, 1.0f, 0.0f);
    glVertex3f(0.25f, 0.75f, 0.0f);
    glVertex3f(-0.25f, 0.75f, 0.0f);
    glVertex3f(0.0f, 1.0f, 0.0f);
    glEnd();
    glBegin(GL_LINE_STRIP);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 1.0f);
    glVertex3f(0.25f, 0.0f, 0.75f);
    glVertex3f(-0.25f, 0.0f, 0.75f);
    glVertex3f(0.0f, 0.0f, 1.0f);
    glVertex3f(0.0f, 0.25f, 0.75f);
    glVertex3f(0.0f, -0.25f, 0.75f);
    glVertex3f(0.0f, 0.0f, 1.0f);
    glEnd();

    glColor4ub(255, 255, 0, 255);
    glRasterPos3f(1.1f, 0.0f, 0.0f);

    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, 'x');
    glRasterPos3f(0.0f, 1.1f, 0.0f);
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, 'y');
    glRasterPos3f(0.0f, 0.0f, 1.1f);
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, 'z');

    glEndList();
}

void display(void)
{
    int i, j, k;
    bool state;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glPushMatrix();

       glTranslatef(0, 0, model_view.z_trans);

       glRotatef(model_view.x_rot, 1, 0, 0);
       glRotatef(model_view.y_rot, 0, 1, 0);

        /* Draw axes */
        glPushMatrix();
            glCallList(axes_list);
        glPopMatrix();


        // Save the lighting state variables
        glPushAttrib(GL_LIGHTING_BIT);
        glDisable(GL_LIGHTING);
        glPushMatrix();
            glColor3f(0,1,0);

            glScalef(WX, WY, WZ);
            glutWireCube(1.0);
        glPopMatrix();
        // Restore lighting state variables
        glPopAttrib();

        for(const auto& p1 : engine->particles) {
            auto p = *p1;
            glPushMatrix();
                glTranslatef(p.pos.x,p.pos.y,p.pos.z);
                SPHEngine::triple h =  engine->phash(p);
                double m = 20.0d;

                if((h.second.second+3)%3 ==1)
                    glColor3f(1,0,0);
                else if((h.second.second+3)%3 ==2)
                    glColor3f(0,1,0);
                else
                    glColor3f(0,0,1);

                //glColor3f(0,0,(50+(double)30*(h.second.second+1))/255.0d);

                //}
               // glutSolidSphere(PRADIUS/2,20,20);
                 glutWireSphere(PRADIUS/2,20,20);
            glPopMatrix();
        }

    glPopMatrix();
    glutSwapBuffers();
}


void simulationRun(void)
{
    //engine->evolve(COLS);
    engine->step(DT);
    glutPostRedisplay();

}

void init(void)
{
    GLfloat  ambientLight[] = { 0.2f, 0.2f, 0.2f, 1.0f };
    GLfloat  diffuseLight[] = { 0.75f, 0.75f, 0.75f, 1.0f };

    glEnable(GL_LIGHTING);
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambientLight);
    glLightfv(GL_LIGHT0, GL_DIFFUSE,  diffuseLight);
    glEnable(GL_LIGHT0);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

    glClearColor (0.0, 0.0, 0.0, 0.0);
    glShadeModel (GL_FLAT);

    glEnable (GL_DEPTH_TEST);

    model_view.x_rot = 0.0;
    model_view.y_rot = 0.0;
    model_view.z_trans = 0.0;


    printf("The life 3D cellular automata model\n");
    printf("Left click on the graphic window to start the simulation\n");
    printf("Right click on the graphic window to stop the simulation\n");
}

void reshape(int w, int h)
{
    GLfloat	 lightPos[]	= { 0.0f, 0.0f, 100.0f, 1.0f };
    double MAX = 4*WX;

    if (MAX < COLS)
        MAX = COLS;
    if (MAX < LAYERS)
        MAX = LAYERS;

    glViewport (0, 0, (GLsizei) w, (GLsizei) h);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    gluPerspective(45.0, (GLfloat) w/(GLfloat) h, 1.0, 8*MAX);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt (0.0, 0.0, 2*MAX, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);

    lightPos[2] = 2*MAX;
    glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
}

void mouse(int button, int state, int x, int y)
{
    switch (button) {
    case GLUT_LEFT_BUTTON:
        if (state == GLUT_DOWN)
        {
            start_time = time(NULL);
            glutIdleFunc(simulationRun);
            //engine->step(DT);
            glutPostRedisplay();

        }
        break;
    case GLUT_MIDDLE_BUTTON:
    case GLUT_RIGHT_BUTTON:
        if (state == GLUT_DOWN)
            glutIdleFunc(NULL);
        break;
    default:
        break;
    }
}

void specialKeys(int key, int x, int y){

    GLubyte specialKey = glutGetModifiers();
    const GLfloat x_rot = 5.0, y_rot = 5.0, z_trans = 0.2;

    if(key==GLUT_KEY_DOWN){
        model_view.x_rot += x_rot;
    }
    if(key==GLUT_KEY_UP){
        model_view.x_rot -= x_rot;
    }
    if(key==GLUT_KEY_LEFT){
        model_view.y_rot -= y_rot;
    }
    if(key==GLUT_KEY_RIGHT){
        model_view.y_rot += y_rot;
    }
    if(key == GLUT_KEY_PAGE_UP){
        model_view.z_trans += z_trans;
    }
    if(key == GLUT_KEY_PAGE_DOWN){
        model_view.z_trans -= z_trans;
    }

    glutPostRedisplay();
}

#define ROWS (0.2d)
#define COLS (0.2d)
#define LAYERS (0.2d)

int main(int argc, char** argv)
{


    engine = new SPHEngine();


    engine->initialize(WX,WY,WZ,MASS,PRADIUS);

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(640, 480);
    glutInitWindowPosition(100, 100);
    glutCreateWindow(argv[0]);
    init();
    initAxesList();
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutSpecialFunc(specialKeys);
    glutMouseFunc(mouse);
    glutMainLoop();

    delete engine;
    std::cout<<"hello"<<std::endl;
    return 0;
}
