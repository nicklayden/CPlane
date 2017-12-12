/*
    Main source file for a c++ plotting interface. 
    Written with SFML.

*/

#include <SFML/Graphics.hpp>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <vector>
#include <string>
#include "plotter.hpp"


Plot::Plot(std::string title, float width, float height, int sizex, int sizey)
:title(title), mainwindow(sf::VideoMode(sizex,sizey),title)
{
    // initialize parameters of the main drawing window.
    plotView.setCenter(0,0);
    plotView.setSize(width,-height);
    plotView.setViewport(sf::FloatRect(0.1,0.1,0.8,0.8));

    axesView.setCenter(0,0);
    axesView.setSize(800,800);
    axesView.setViewport(sf::FloatRect(0,0,1,1));
    mainwindow.clear(sf::Color::Black);

    // Load font package
    // font.loadFromFile("OpenSans-Regular.ttf");


    // Initialize text regions.
    // texttitle.setFont(font);
    // texttitle.setString("Test");
    // texttitle.setCharacterSize(12);
    // texttitle.setColor(sf::Color::White);
    // texttitle.setPosition(0,0);

}

void Plot::plot(std::vector<double> x, std::vector<double> y, sf::Color colour)
{
    // Check vector bounds to rescale axes. only if autoscale = true
    // CheckBounds(x,y);
    
    EventLoop();
    DrawLines(x,y, colour);
    // mainwindow.display();
}

void Plot::scatter(std::vector<double> x, std::vector<double> y, sf::Color colour)
{
    EventLoop();
    DrawPoints(x,y,colour);

}

void Plot::scatter(double x, double y, sf::Color colour)
{
    EventLoop();
    sf::CircleShape pt(0.01);
    pt.setOrigin(pt.getRadius(),pt.getRadius());
    pt.setPosition(x,y);
    pt.setFillColor(colour);
    mainwindow.draw(pt);
}

void Plot::show()
{
    mainwindow.display();
    mainwindow.clear(sf::Color::Black);
}

void Plot::DrawLines(std::vector<double> x, std::vector<double> y, sf::Color colour)
{
    int xsize = x.size(); int ysize = y.size();
    if (xsize == ysize)
    {
        // Draws lines between a series of points.
        sf::VertexArray solidline(sf::LinesStrip,x.size());
            for (size_t i = 0; i < xsize; i++) {
                solidline[i].position = sf::Vector2f(x[i],y[i]);
                // if (i%2 == 0) {
                solidline[i].color = colour;
                // } else {
                solidline[i].color = colour;
                // }
            }
            mainwindow.draw(solidline);
    } else {
        std::cerr << "ERROR: x and y must be the same size.\n";
    }
}

void Plot::DrawPoints(std::vector<double> x, std::vector<double> y, sf::Color colour)
{
    sf::CircleShape pt(0.01);
    for (size_t i = 0; i < x.size(); i++) {
        pt.setOrigin(pt.getRadius(), pt.getRadius());
        pt.setPosition(x[i],y[i]);
        pt.setFillColor(colour);
        mainwindow.draw(pt);
    }

}

void Plot::EventLoop()
{
    sf::Event event;
    while (mainwindow.pollEvent(event))
    {
        if (event.type == sf::Event::Closed)
        {
            mainwindow.close();
        }   
    }
    // mainwindow.clear(sf::Color::Black);
    // Draw axes bounds and (eventually) labels
    mainwindow.setView(axesView);
    DrawBoundingBox(sf::Color::White);
    // mainwindow.draw(texttitle);

    mainwindow.setView(plotView);
    // double width = abs(xmin) + abs(xmax);
    // double height= abs(ymin) + abs(ymax);
    // plotView.setSize(width,-height);
    // mainwindow.display();
}

void Plot::DrawBoundingBox(sf::Color colour)
{
    // Draws the bounding box for the axes around the window.
    sf::VertexArray boundingbox(sf::LinesStrip, 5);
    xmin = -320; xmax = 320;
    ymin = -320; ymax = 320;
    boundingbox[0].position = sf::Vector2f(xmin, ymax);
    boundingbox[1].position = sf::Vector2f(xmax, ymax);
    boundingbox[2].position = sf::Vector2f(xmax, ymin);
    boundingbox[3].position = sf::Vector2f(xmin, ymin);
    boundingbox[4].position = sf::Vector2f(xmin, ymax);

    for (size_t i = 0; i < 5; i++) {
        boundingbox[i].color = colour;
    }
    mainwindow.draw(boundingbox);
}

inline void Plot::CheckBounds(std::vector<double> x, std::vector<double> y)
{
    auto xbounds = std::minmax_element(x.begin(),x.end());
    auto ybounds = std::minmax_element(y.begin(),y.end());
    std::cout << "X bounds: [" << *xbounds.first << "," << *xbounds.second << "]" << std::endl;
    std::cout << "Y bounds: [" << *ybounds.first << "," << *ybounds.second << "]" << std::endl;
}

void Plot::setxlim(double xmin, double xmax)
{
    this->xmin = xmin;
    this->xmax = xmax;

}

void Plot::setylim(double ymin, double ymax)
{
    this->ymin = ymin;
    this->ymax = ymax;

}

void Plot::setTitle(std::string title)
{

}

std::string Plot::toString(double number)
{

}

void Plot::set_xlabel(std::string label, sf::Font font)
{

    xlabel.setFont(font);
    xlabel.setString(label);
    xlabel.setPosition(-320,-350);
    xlabel.setCharacterSize(16);
    xlabel.setColor(sf::Color::White);
    mainwindow.draw(xlabel);
}

void Plot::set_ylabel(std::string)
{
    
}

bool Plot::isOpen()
{
    return mainwindow.isOpen();
}