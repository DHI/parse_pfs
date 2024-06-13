import mikeio
import pandas as pd
import os
from pyproj import CRS
from mikeio.pfs._pfssection import PfsSection
from datetime import datetime, timedelta

class PFS_Parser:
    def __init__(self, fname):
        self.fname = fname
        self.pfs = mikeio.read_pfs(fname)["FemEngineHD"]
        self.domain = self.pfs["DOMAIN"]
        self.time = self.pfs["TIME"]
        self.hydrodynamic = self.pfs["HYDRODYNAMIC_MODULE"]

    def parse(self):
        self.attributes = {}
        self.parse_M3FM_Domain()
        self.parse_M3FM_Time()
        self.parse_M3FM_Hydrodynamic()
        
        return self.attributes

    def parse_M3FM_Domain(self):
        """
        Parse the domain section of a PFS file
        """
        # Coordinate system
        try:
            crs = CRS.from_string(self.domain["coordinate_type"])
            coordinate_system = crs.name
        except:
            coordinate_system = self.domain["coordinate_type"]
        # Vertical discretization
        if self.domain["vertical_mesh_type_overall"] == 1:
            # σ-coordinate
            number_of_layers = self.domain["number_of_layers"]
            if self.domain["vertical_mesh_type"] == 1:
                mesh_type = "equidistant layers"
            elif self.domain["vertical_mesh_type"] == 2:
                mesh_type = "layers with specified thicknesses"
            elif self.domain["vertical_mesh_type"] == 3:
                mesh_type = "layers with specified distribution parameters"
            vertical_mesh_type = f"σ-coordinate with {number_of_layers} {mesh_type}"
        elif self.domain["vertical_mesh_type_overall"] == 2:
            # σ-coordinate
            number_of_layers = self.domain["number_of_layers"]
            if self.domain["vertical_mesh_type"] == 1:
                mesh_type = "equidistant layers"
            elif self.domain["vertical_mesh_type"] == 2:
                mesh_type = "layers with specified thicknesses"
            elif self.domain["vertical_mesh_type"] == 3:
                mesh_type = "layers with specified distribution parameters"
            top_vertical_mesh = f"σ-coordinate with {number_of_layers} {mesh_type} at the top"
            # z-coordinate
            number_of_layers = self.domain["number_of_layers_zlevel"]
            if self.domain["vertical_mesh_type_zlevel"] == 1:
                mesh_type = "equidistant layers"
            elif self.domain["vertical_mesh_type_zlevel"] == 2:
                mesh_type = "layers with specified thicknesses"
            bottom_vertical_mesh = f"z-coordinate with {number_of_layers} {mesh_type} at the bottom"
            vertical_mesh_type = f"{top_vertical_mesh} and {bottom_vertical_mesh}"
        
        # Boundary names
        boundary_names = self._parse_Boundaries(self.domain["BOUNDARY_NAMES"])
        
        self.attributes = self.attributes | {"Domain":
                                             {"Coordinate System": coordinate_system,
                                              "Vertical Mesh": vertical_mesh_type,
                                              "Boundaries": boundary_names}
                                            }

    def parse_M3FM_Time(self):
        """
        Parse the time section of a PFS file
        """
        start_time = datetime(*self.time["start_time"])
        time_step_interval = f'{self.time["time_step_interval"]} s'
        number_of_time_steps = self.time["number_of_time_steps"]
        end_time = start_time + timedelta(seconds=self.time["number_of_time_steps"]*self.time["time_step_interval"])
        self.attributes = self.attributes | {"Time":
                                             {"Start Time": start_time.strftime("%Y-%m-%d %I:%M:%S %p"),
                                              "Time Step Interval": time_step_interval,
                                              "Number of Time Steps": str(number_of_time_steps),
                                              "End Time": end_time.strftime("%Y-%m-%d %I:%M:%S %p")}
                                            }

    def parse_M3FM_Hydrodynamic(self):
        if self.hydrodynamic["mode"] != 2:
            return "Hydrodynamic mode is not activated!"
        governing_equations = "Navier-Stokes equations" if self.hydrodynamic["EQUATION"]["formulation"] == 2 else "Shallow water equations" if self.hydrodynamic["EQUATION"]["formulation"] == 4 else None
        governing_equations = {"Governing Equation": governing_equations}
        solution_technique = self._parse_Solution_Technique(self.hydrodynamic["SOLUTION_TECHNIQUE"])
        depth = self._parse_Depth(self.hydrodynamic["DEPTH"])
        flood_and_dry = self._parse_Flood_Dry(self.hydrodynamic["FLOOD_AND_DRY"])
        density = self._parse_density(self.hydrodynamic["DENSITY"])
        turbulence = self._parse_Turbulence(self.hydrodynamic["EDDY_VISCOSITY"])
        bed_resistence = self._parse_Bed_Resistence(self.hydrodynamic["BED_RESISTANCE"])
        vegetation = "No vegetation" if self.hydrodynamic["VEGETATION"]["type"] == 0 else "Vegetation Included"
        coriolis = self._parse_Coriolis(self.hydrodynamic["CORIOLIS"])
        wind = self._parse_Wind_Forcing(self.hydrodynamic["WIND_FORCING"])
        ice = self._parse_Ice_Coverage(self.hydrodynamic["ICE"])
        tide = self._parse_Tidal_Potential(self.hydrodynamic["TIDAL_POTENTIAL"])
        self.attributes = self.attributes | {"Hydrodynamics":
                                             {"Numerical Characteristics": governing_equations | solution_technique,
                                              "Depth Correction Type": depth,
                                              "Flood and Dry Type": flood_and_dry,
                                              "Density": density,
                                              "Turbulence": turbulence,
                                              "Bed Resistence": bed_resistence,
                                              "Vegetation": vegetation,
                                              "Coriolis": coriolis,
                                              "Wind Forcing": wind,
                                              "Ice Coverage": ice,
                                              "Tidal Potential": tide}
                                            }

    def _parse_Boundaries(self, boundaries):
        self.boundaries = []
        for key in boundaries.keys():
            if isinstance(boundaries[key], PfsSection):
                self.boundaries.append(boundaries[key]["name"])
        
    def _parse_Solution_Technique(self, solution_technique):
        schemes = {1: "Low order, fast algorithm", 2: "Higher order", 3: "Higher order"}
        Reimann_solvers = {0: None, 3: "HLLC"}
        Filtering_options = {0: None, 1: "Explicit filter"}
        # HD equations
        time_integration_scheme = schemes[solution_technique["scheme_of_time_integration"]]
        horizontal_discretization = schemes[solution_technique["scheme_of_space_discretization_horizontal"]]
        vertical_discretization = schemes[solution_technique["scheme_of_space_discretization_vertical"]]
        dt_max = solution_technique["dt_max_HD"]
        dt_min = solution_technique["dt_min_HD"]
        critical_CFL = solution_technique["CFL_critical_HD"]
        HD = {"Hydrodynamic Solver": {"Time Integration Scheme": time_integration_scheme, "Horizontal Discretization": horizontal_discretization, "Vertical Discretization": vertical_discretization, "Maximum Time Step": dt_max, "Minimum Time Step": dt_min, "Critical CFL": critical_CFL}}
        # Sediment transport equations
        dt_max = solution_technique["dt_max_AD"]
        dt_min = solution_technique["dt_min_AD"]
        critical_CFL = solution_technique["CFL_critical_AD"]
        AD = {"Sediment Transport Solver": {"Maximum Time Step": dt_max, "Minimum Time Step": dt_min, "Critical CFL": critical_CFL}}
        # Reimann solver
        Reimann_solver = Reimann_solvers[solution_technique["type_of_Riemann_solver"]]
        Reimann_factor = solution_technique["Riemann_factor"] if Reimann_solver is not None else None
        Reimann = {"Reimann": {"Reimann Solver": Reimann_solver, "Reimann Factor": Reimann_factor}}
        # Filteration
        Filter_type = Filtering_options[solution_technique["type_of_filtering"]]
        Filter_coefficient = solution_technique["filtering_coefficient"] if Filter_type is not None else None
        Filter = {"Filtering": {"Filtering Type": Filter_type, "Filtering Coefficient": Filter_coefficient}}
        return HD | AD | Reimann | Filter

    def _parse_Depth(self, depth):
        # Depth
        if depth["type"] == 0:
            return "No depth correction"
        else:
            if depth["format"] == 2:
                return "Varying in space, constant in time depth correction"
            elif depth["format"] == 3:
                return "Varying in time and domain depth correction"

    def _parse_Flood_Dry(self, flood_and_dry):
        if flood_and_dry["type"] == 0:
            return {}
        elif flood_and_dry["type"] == 2:
            dry = flood_and_dry["drying_depth"]
            wet = flood_and_dry["mass_depth"]
            return {"Drying depth": dry, "Wetting depth": wet}

    def _parse_density(self, density):
        if density["type"] == 0:
            return {"Density Type": "Barotropic"}
        elif density["type"] == 1:
            return {"Density type": "Function of temperature and salinity", "Reference Temperature": density["temperature_reference"], "Reference Salinity": density["salinity_reference"]}
        elif density["type"] == 2:
            return {"Density type": "Function of temperature", "Reference Temperature": density["temperature_reference"]}
        elif density["type"] == 3:
            return {"Density type": "Function of salinity", "Reference Salinity": density["salinity_reference"]}

    def _parse_Turbulence(self, turbulence):
        horizontal = turbulence["HORIZONTAL_EDDY_VISCOSITY"]
        vertical = turbulence["VERTICAL_EDDY_VISCOSITY"]
        functions = {1: self._parse_Constant_Eddy, 3: self._parse_Smagorinsky_Eddy, 4: self._parse_Log_Law_Eddy, 5: self._parse_Two_Equation_Eddy}
        if horizontal["type"] == 0:
            horizontal_eddy = {"Horizontal Eddy Viscosity": "No eddy"}
        else:
            horizontal_eddy = {"Horizontal Eddy Viscosity": functions[horizontal["type"]](horizontal)}
        if vertical["type"] == 0:
            vertical_eddy = {"Vertical Eddy Viscosity": "No eddy"}
        else:
            vertical_eddy = {"Vertical Eddy Viscosity": functions[vertical["type"]](vertical)}
        return horizontal_eddy | vertical_eddy

    def _parse_Constant_Eddy(self, eddy):
        eddy = eddy["CONSTANT_EDDY_FORMULATION"]
        if eddy["format"] == 0:
            part1 = {"Eddy Viscosity Type": "Constant eddy formulation", "Eddy Viscosity": eddy["constant_value"]}
        elif eddy["format"] == 2:
            part1 = {"Eddy Viscosity Type": "Varying in space"}
        try:
            if eddy["Ri_damping"] == 1:
                part2 = {"Damping constant a": eddy["Ri_a"], "Damping constant b": eddy["Ri_b"]}
        except:
            part2 = {}
        return part1 | part2

    def _parse_Smagorinsky_Eddy(self, eddy):
        eddy = eddy["SMAGORINSKY_FORMULATION"]
        if eddy["format"] == 0:
            part1 = {"Eddy Viscosity Type": "Smagorinsky formulation", "Smagorinsky constant": eddy["constant_value"]}
        elif eddy["format"] == 2:
            part1 = {"Eddy Viscosity Type": "Varying in space"}
        part2 = {"Minimum eddy viscosity": eddy["minimum_eddy_viscosity"], "Maximum eddy viscosity": eddy["maximum_eddy_viscosity"]}
        return part1 | part2
    
    def _parse_Two_Equation_Eddy(self, eddy):
        eddy = eddy["K_EPSILON_FORMULATION"]
        part1 = {"Eddy Viscosity Type": "Two-equation turbulence formulation", "Minimum eddy viscosity": eddy["minimum_eddy_viscosity"], "Maximum eddy viscosity": eddy["maximum_eddy_viscosity"]}
        return part1

    def _parse_Log_Law_Eddy(self, eddy):
        eddy = eddy["LOG_LAW_FORMULATION"]
        part1 = {"Eddy Viscosity Type": "Log-law formulation", "Minimum eddy viscosity": eddy["minimum_eddy_viscosity"], "Maximum eddy viscosity": eddy["maximum_eddy_viscosity"]}
        try:
            if eddy["Ri_damping"] == 1:
                part2 = {"Damping constant a": eddy["Ri_a"], "Damping constant b": eddy["Ri_b"]}
        except:
            part2 = {}
        return part1 | part2

    def _parse_Bed_Resistence(self, bed_resistence):
        if bed_resistence["type"] == 0:
            return {"Bed Resistence Type": "No bed resistence"}
        elif bed_resistence["type"] == 2:
            return {"Bed Resistence Type": "Quadratic drag coefficient"} | self._parse_Quadratic_Drag_Resistence(bed_resistence["DRAG_COEFFICIENT"])
        elif bed_resistence["type"] == 5:
            return {"Bed Resistence Type": "Roughness height"} | self._parse_Roughness_Height_Resistence(bed_resistence["ROUGHNESS"])
        elif bed_resistence["type"] == 7:
            return {"Bed Resistence Type": "Wave induced bed resistence"} | self._parse_Wave_Induced_Resistence(bed_resistence["WAVE_INDUCED_ROUGHNESS"])

    def _parse_Quadratic_Drag_Resistence(self, drag_coefficient):
        if drag_coefficient["format"] == 0:
            return {"Bed resistence format": "Constant", "Constant value": drag_coefficient["constant_value"]}
        elif drag_coefficient["format"] == 2:
            return {"Bed resistence format": "Constant in time, varying in space"}
        elif drag_coefficient["format"] == 3:
            return {"Bed resistence format": "Varying in time and space"}
        
    def _parse_Roughness_Height_Resistence(self, roughness_height):
        if roughness_height["format"] == 0:
            return {"Bed resistence format": "Constant", "Constant value": roughness_height["constant_value"]}
        elif roughness_height["format"] == 2:
            return {"Bed resistence format": "Constant in time, varying in space"}
        elif roughness_height["format"] == 3:
            return {"Bed resistence format": "Varying in time and space"}
        
    def _parse_Wave_Induced_Resistence(self, wave_data):
        if wave_data["format"] == 0:
            part1 = {"Bed resistence format": "Constant", "Constant value": wave_data["constant_value"]}
        elif wave_data["format"] == 2:
            part1 = {"Bed resistence format": "Constant in time, varying in space"}
        elif wave_data["format"] == 3:
            part1 = {"Bed resistence format": "Varying in time and space"}

        if wave_data["type_of_wave_height_correction"] == 0:
            part2 = {"Wave correction type": "No correction"}
        elif wave_data["type_of_wave_height_correction"] == 1:
            part2 = {"Wave correction type": "Wave height correction", "Wave height over depth limit": wave_data["wave_height_over_depth_limit"], "Minimum water depth for including waves": wave_data["minimum_waterdepth_for_including_waves"]}

        return part1 | part2

    def _parse_Coriolis(self, coriolis):
        if coriolis["type"] == 0:
            return {"Coriolis forcing": "No coriolis force"}
        elif coriolis["type"] == 1:
            return {"Coriolis forcing": "Constant coriolis forcing in space", "Reference latitude": coriolis["latitude"]}
        elif coriolis["type"] == 2:
            return {"Coriolis forcing": "Varying coriolis forcing in space"}

    def _parse_Wind_Forcing(self, wind):
        if wind["type"] == 0:
            return {"Wind forcing type": "No wind forcing"}
        elif wind["type"] == 1:
            if wind["format"] == 0:
                part1 = {"Wind forcing format": "Constant", "Constant wind speed": wind["constant_speed"], "Constant wind direction": wind["constant_direction"], "Soft time interval": str(wind["soft_time_interval"]) + " s"}
            elif wind["format"] == 1:
                part1 = {"Wind forcing format": "Varying in time, constant in space", "Soft time interval": str(wind["soft_time_interval"]) + " s"}
            elif wind["format"] == 3:
                part1 = {"Wind forcing format": "Varying in time and space", "Neutral pressure": str(wind["neutral_pressure"]) + " hPa", "Soft time interval": str(wind["soft_time_interval"]) + " s"}
            wind_friction = wind["WIND_FRICTION"]
            if wind_friction["type"] == 0:
                part2 = {"Wind friction type": "Constant", "Constant wind friction": wind_friction["constant_friction"]}
            elif wind_friction["type"] == 1:
                lower_speed = wind_friction["linear_speed_low"]; upper_speed = wind_friction["linear_speed_high"]
                lower_friction = wind_friction["linear_friction_low"]; upper_friction = wind_friction["linear_friction_high"]
                part2 = {"Wind friction type": "Function of wind speed", "Lower limit for linear variation": f"{lower_friction} at {lower_speed} m/s", "Upper limit for linear variation": f"{upper_friction} at {upper_speed} m/s"}
            return part1 | part2

    def _parse_Ice_Coverage(self, ice):
        if ice["type"] == 0:
            return {"Ice coverage type": "No ice coverage"}
        else:
            if ice["type"] == 1:
                part1 = {"Ice coverage type": "Specified ice concentration", "Critical ice concentration": ice["c_cut_off"]}
            elif ice["type"] == 2:
                part1 = {"Ice coverage type": "Specified ice thickness"}
            elif ice["type"] == 3:
                part1 = {"Ice coverage type": "Specified ice concentration and thickness", "Critical ice concentration": ice["c_cut_off"]}
            ice_roughness = ice["ROUGHNESS"]
            if ice_roughness["type"] == 0:
                part2 = {"Ice roughness": "Not included"}
            else:
                if ice_roughness["format"] == 0:
                    part2 = {"Ice roughness type": "Constant", "Constant ice roughness height": str(ice_roughness["constant_value"]) + " m"}
                elif ice_roughness["format"] == 2:
                    part2 = {"Ice roughness type": "Varying in space"}
            return part1 | part2

    def _parse_Tidal_Potential(self, tide):
        if tide["type"] == 0:
            return {"Tidal potential": "Not included"}
        else:
            if tide["format"] == 0:
                n_constituents = tide["number_of_constituents"]
                constituents = " ".join([tide[f"CONSTITUENT_{i+1}"]["name"] for i in range(n_constituents)])
                return {"Tidal potential": "Included", "Number of constituents": n_constituents, "Constituents": constituents}
            elif tide["format"] == 1:
                return {"Tidal potential": "Included"}

sample_files = [f for f in os.listdir("sample_input_files")]
fname = os.path.join("sample_input_files", sample_files[9])

pfs = PFS_Parser(fname)
attributes = pfs.parse()
print(attributes)