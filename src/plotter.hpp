#pragma once
#include <SFML/Graphics.hpp>
#include <vector>
#include <string>
#include <algorithm>


/*
    Input values to plot as vectors. fuck off with the other data types.
*/
class Plot
{
    public:
        Plot(std::string title, float width = 5, float height = 5, int sizex = 800, int sizey = 800);
        void plot(std::vector<double> x, std::vector<double> y,  sf::Color color = sf::Color::Blue); 
        void scatter(std::vector<double> x, std::vector<double> y,  sf::Color color = sf::Color::Blue);
        void scatter(double x, double y, sf::Color colour = sf::Color::Blue);
        void show();

        sf::RenderWindow mainwindow;
        // sf::Text texttitle;
        double xmax = .5;
        double xmin = -.5;
        double ymin = -.5;
        double ymax = .5;

        void EventLoop();
        void animate();
        void setTitle(std::string title);
        std::string toString(double number);
        void setxlim(double xmin, double xmax);
        void setylim(double ymin, double ymax);
        void set_xlabel(std::string label, sf::Font font);
        void set_ylabel(std::string);
        sf::View axesView;
        sf::View plotView;

        sf::Text xlabel;
        sf::Text ylabel;

    private:
        void DrawLines(std::vector<double> x, std::vector<double> y, sf::Color colour);
        void DrawPoints(std::vector<double> x, std::vector<double> y, sf::Color colour);
        void DrawBoundingBox(sf::Color colour);
        void CheckBounds(std::vector<double> x, std::vector<double> y);

        std::string title;
        bool autoscale = false;
        float xscale = 1.2;
        float yscale = 1.2;
};