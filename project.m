% Merry-Go-Round Project
% Matthew Krawczyk
% Karl Auer
%
% AE138 Vector Dynamics
% 6 December 2019
%------------------------------------------------------------------------------------------------------------------
% Author Note: this code was written to be run as a MATLAB(TM) script. Do not attempt to use this 
%              code in MotionGenesis
%------------------------------------------------------------------------------------------------------------------

% Problem statement:
% A rider is on a simple merry-go-round at some amusement park
% The ride consists of an outer platform of radius R and a smaller "cup" (with radius r, and distance d away from 
% the center of the platform) which houses the riders
% 
% Part 1:
%   Find the amount of time it would take a rider to reach a perceived acceleration of 6g (~58.8 m/s) if both
%   platforms start at rest and are continuously accelerated at a constant rate.
%
% Part 2:
%   Determine the angular velocity at which both the platform and the "cup" would be spinning at when the rider 
%   hits a perceived acceleration of 6g.
%
% Assumptions:
%   -The distance from the center of the platform to the center of the cup does not change
%   -The rider stays a constant distance r from the center of the cup
%   -The effects of air resistance on the rider are neglected
%   -The rider can be treated as a point mass
%   
%   -Acceleration of the platform and the cup is constant
%   -Main platform and "cup" start at rest
%   -The main platform and rider begin oriented horizontally

clc;                                    % clear screen
close all;                              % close any open figures
%------------------------------------------------------------------------------------------------------------------
% Define the constants
distance = 4;                           % Distance from the center of the platform to the center of the cup (in m)
radius = 0.5;                           % Distance from the center of the cup to the rider (in m)
theta_dPrime = 0.1;                     % Angular acceleration of the main platfrom (in deg/s^2)
phi_dPrime = 0.3;                       % Angular acceleration of the "cup" (in deg/s^2)
gravity = 9.8;                          % Gravity (9.8 m/s^2)
max_g = 6;                              % Maximum g-force before blackout (for an average person)

step_size = 0.01;                       % Time step for computations

%------------------------------------------------------------------------------------------------------------------
% Call the functions and print out result
max_t = part1(distance, radius, theta_dPrime, phi_dPrime, gravity, max_g, step_size);
fprintf('It will take %.2fs for the rider to experience %.0fg''s.\n', max_t, max_g);
fprintf('The main platform will be spinning at %.2f deg/s and the cup will be spinning at %.2f deg/s.\n', ...
    theta_dPrime * max_t, phi_dPrime * max_t);
graphResult(distance, radius, theta_dPrime, phi_dPrime, max_t, gravity, max_g, step_size);

%------------------------------------------------------------------------------------------------------------------
% Part 1
% How long to hit 6g's?
% @param distance - distance from the center of the platform to the center of the cup
% @param radius - distance from the center of the cup to the rider
% @param theta_dPrime - angular acceleration of the main platform
% @param phi_dPrime - angular acceleration of the cup
% @param gravity - force of gravity
% @param max_g - maximum g-force before blacking out
% @param step_size - the granularity of the data (how large each time step is)
%
% @return max_time - time computed to hit {@param max_g}
function max_time = part1(distance, radius, theta_dPrime, phi_dPrime, gravity, max_g, step_size)
    increment = 0;                      % Counter for the loop
    
    while (1)
        t = increment * step_size;      % time
    
        % NB: with a constant acceleration x (and starting with no angular velocity and oriented horizontally)
        % the velocity (w.r.t time) is x*t and the position (w.r.t time) is 0.5*x*t^2
        phi = 0.5 * phi_dPrime * t^2;
        theta_prime = theta_dPrime * t;
        phi_prime = phi_dPrime * t;
    
        % x and y components of acceleration
        accel_x = (-(distance * theta_dPrime) + (sind(phi) * (theta_dPrime + phi_dPrime) * radius) - ...
            (((theta_prime + phi_prime)^2) * radius * cosd(phi)));
        accel_y = (-(distance * theta_prime^2) - (cosd(phi) * (theta_dPrime + phi_dPrime) * radius) - ...
            (((theta_prime + phi_prime)^2) * radius * sind(phi)));
    
        % total acceleration
        accel = sqrt(accel_x^2 + accel_y^2);
    
        increment = increment + 1;
    
        if (accel >= (max_g * gravity))
            max_time = t;
            return;
        else
            increment = increment + 1;
        end
    end
end

%------------------------------------------------------------------------------------------------------------------
% Show the results in a graphical format
% @param distance - distance from the center of the platform to the center of the cup
% @param radius - distance from the center of the cup to the rider
% @param theta_dPrime - angular acceleration of the main platform
% @param phi_dPrime - angular acceleration of the cup
% @param max_t - time computed to hit {@param max_g}, returned by part1()
% @param gravity - force of gravity
% @param max_g - maximum g-force before blacking out
% @param step_size - the granularity of the data (how large each time step is)
%
% @note - {@param max_g} and {@param gravity} are used to draw a referece line in the acceleration
%
% @print figure1 - A graph of the displacement vs. time
% @print figure2 - A graph of the angular velocity vs. time
% @print figure3 - A graph of the acceleration vs. time
function graphResult(distance, radius, theta_dPrime, phi_dPrime, max_t, gravity, max_g, step_size)
    % Preallocate array sizes to speed up code
    time = zeros(1, (max_t / step_size) + 1);
    % Angles
    theta_t = zeros(1, (max_t / step_size) + 1);
    phi_t = zeros(1, (max_t / step_size) + 1);
    % Accelerations
    accel_x_t = zeros(1, (max_t / step_size) + 1);
    accel_y_t = zeros(1, (max_t / step_size) + 1);
    accel_t = zeros(1, (max_t / step_size) + 1);
    
    max_g_t = ones(1, (max_t / step_size) + 1) * (gravity * max_g);
    % Distance from center
    distance_t = zeros(1, (max_t / step_size) + 1);
    
    increment = 0;                      % Counter for the loop
    
    while(1)
        t = increment * step_size;      % time
               
        % NB: with a constant acceleration x (and starting with no angular velocity and oriented horizontally)
        % the velocity (w.r.t time) is x*t and the position (w.r.t time) is 0.5*x*t^2
        phi = 0.5 * phi_dPrime * t^2;
        theta_prime = theta_dPrime * t;
        phi_prime = phi_dPrime * t;
        
        % x and y componenents of acceleration
        accel_x = (-(distance * theta_dPrime) + (sind(phi) * (theta_dPrime + phi_dPrime) * radius) - ...
            (((theta_prime + phi_prime)^2) * radius * cosd(phi)));
        accel_y = (-(distance * theta_prime^2) - (cosd(phi) * (theta_dPrime + phi_dPrime) * radius) - ...
            (((theta_prime + phi_prime)^2) * radius * sind(phi)));
        
        % Distance from center
        distance_x = distance + (radius * cosd(phi));
        distance_y = (radius * sind(phi));
        
        % Save data to arrays
        %time
        time(increment + 1) = t;
        % angles
        theta_t(increment + 1) = theta_dPrime * t;
        phi_t(increment + 1) = phi_dPrime * t;
        % distance
        distance_t(increment + 1) = sqrt(distance_x^2 + distance_y^2);
        % acceleration
        accel_x_t(increment + 1) = accel_x;
        accel_y_t(increment + 1) = accel_y;
        accel_t(increment + 1) = sqrt(accel_x^2 + accel_y^2);
        
        if (t >= max_t)
            break
        else
            increment = increment + 1;
        end
    end
    
    % Displacement curve
    figure();
    plot(time, distance_t, 'b.');
    title('Displacement v. time');
    xlabel('Time (s)');
    ylabel('Displacement (m)');
    legend('Distance from the center to the rider');
    
    % Velocity curve
    figure();
    plot(time, theta_t, 'k.', time, phi_t, 'b.');
    title('Angular velocty v. time');
    xlabel('Time (s)');
    ylabel('Angular velocity (deg/s)');
    legend('Main Platform', 'Cup');
    
    % Acceleration curve
    figure();
    plot(time, accel_t, 'k.', time, accel_x_t, 'b.', time, accel_y_t, 'g.', time, max_g_t, 'r.');
    title('Acceleration v. Time');
    xlabel('Time (s)');
    ylabel('Acceleration (m/s^2)');
    legend('Magnitude of total acceleration', 'x-component of acceleration', 'y-component of acceleration', 'Maximum acceleration before g-LOC');
end
