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
        Plot(std::string title);
        void plot(std::vector<double> x, std::vector<double> y,  sf::Color color = sf::Color::Blue); 
        void show();

        sf::RenderWindow mainwindow;
        // sf::Text title;

        double xmax = 15;
        double xmin = -15;
        double ymin = -15;
        double ymax = 15;

        void EventLoop();
        void animate();
    private:
        void DrawLines(std::vector<double> x, std::vector<double> y, sf::Color colour);
        void DrawPoints(std::vector<double> x, std::vector<double> y);
        void DrawBoundingBox(sf::Color colour);
        void CheckBounds(std::vector<double> x, std::vector<double> y);

        sf::View plotView;
        sf::View axesView;
        std::string title;
        bool autoscale = false;
        float xscale = 1.2;
        float yscale = 1.2;
};