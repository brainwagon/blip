I need assistance in developing a Python program to analyze and optimize the placement of three symmetric support points for a mirror cell in an amateur Newtonian telescope. The program should use finite element analysis (FEA) to evaluate the deformation of the primary mirror under gravitational load and determine the root mean square (RMS) deviation of the mirror surface. Additionally, it should identify the optimal support radius that minimizes the RMS deviation. Below are the detailed requirements and specifications for the program:

1. Input Parameters
Mirror Dimensions:
Diameter of the primary mirror (in millimeters).
Thickness of the primary mirror (in millimeters).
Material Properties:
Assume the mirror material is Pyrex (borosilicate glass) with the following properties:
Young's Modulus: 62.7 GPa (or 6.27e10 Pa).
Poisson's Ratio: 0.2.
Density: 2230 kg/m³.
Support Configuration:
Three support points, symmetrically placed at 120-degree intervals around the mirror's center.
Support radius (distance from the center to each support point) as a variable to be optimized or user-specified for verification.
Loading Condition:
Gravitational load acting perpendicular to the mirror surface (assume the telescope is pointing at zenith, so gravity acts along the optical axis).
Acceleration due to gravity: 9.81 m/s².
2. Program Objectives
Calculate Deformation:
Use finite element analysis to model the mirror as a circular plate and compute the deformation (displacement) of the mirror surface under gravitational load for a given support radius.
Compute RMS Deviation:
Calculate the root mean square (RMS) deviation of the mirror surface from its ideal shape (flat or parabolic, approximated as flat for small deformations) at a given support radius.
Express RMS deviation in nanometers or micrometers for optical relevance.
Optimize Support Radius:
Iteratively test different support radii (e.g., from 0.2 * radius to 0.8 * radius) to find the radius that minimizes the RMS deviation.
Output the optimal support radius and the corresponding RMS deviation.
3. Technical Requirements
Finite Element Analysis Library:
Use a Python-compatible FEA library such as FEniCS, scikit-fem, or meshio combined with a solver like Gmsh for mesh generation.
If a full FEA implementation is complex, consider a simplified plate-bending model based on Kirchhoff-Love plate theory as a fallback, with appropriate references or validation.
Mesh Generation:
Generate a 2D or 3D mesh of the circular mirror with sufficient resolution to capture deformation accurately (e.g., triangular or quadrilateral elements).
Account for symmetry to reduce computational load if possible.
Boundary Conditions:
Model the three support points as fixed or pinned constraints (no vertical displacement at support locations).
No additional edge support (free edges) unless specified by the user.
Output Visualization:
Plot the deformation profile of the mirror surface (e.g., using matplotlib or plotly for 2D/3D visualization).
Display the RMS deviation as a function of support radius to show the optimization process.
4. Deliverables
Python Script:
A well-documented Python script with clear comments explaining each section (input, FEA setup, deformation calculation, RMS computation, optimization, and output).
Include error handling for invalid inputs (e.g., negative dimensions, unrealistic support radii).
Sample Output:
Run the program for a sample mirror (e.g., 150mm diameter, 25mm thickness) and provide:
RMS deviation at a user-specified support radius.
Optimal support radius and corresponding minimum RMS deviation.
Visual plot of deformation and/or RMS deviation vs. support radius.
Dependencies:
List all required Python libraries (e.g., numpy, scipy, matplotlib, FEA-specific libraries) and provide installation instructions.
5. Additional Notes
Accuracy and Validation:
Validate the program against known results or software like PLOP (a popular mirror cell optimization tool) if possible.
Ensure the deformation results are physically realistic (e.g., compare RMS deviation to typical optical tolerances of lambda/10 or lambda/20 for visible light, where lambda = 550nm).
Scalability:
Design the code to handle mirrors of different sizes and thicknesses without significant modification.
User-Friendly Interface:
Allow users to input parameters via command line or a simple GUI (e.g., using tkinter) if feasible.
Performance:
Optimize the code for reasonable computation time on a standard desktop computer, balancing mesh resolution with speed.
6. Context and Motivation
This program is intended for amateur telescope makers who need to design mirror cells with minimal optical deformation. The goal is to provide a free, open-source tool that replicates the functionality of existing software like PLOP, focusing on ease of use and transparency in calculations. The program should empower users to verify support placement and understand the impact of support radius on mirror performance.

Please provide the complete Python code, along with explanations of key sections, installation instructions for dependencies, and a sample run for a typical mirror. If full FEA is not feasible within the scope, suggest a simplified approach with references to plate theory or existing libraries that can be adapted. Thank you! 🪞
