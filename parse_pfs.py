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
        coordinate_system = self.__parse_coordinates(self.domain["coordinate_type"])
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
        precipitation = self._parse_Precipitation(self.hydrodynamic["PRECIPITATION_EVAPORATION"]["PRECIPITATION"])
        evaporation = self._parse_Evaporation(self.hydrodynamic["PRECIPITATION_EVAPORATION"]["EVAPORATION"])
        infiltration = self._parse_Infiltration(self.hydrodynamic["INFILTRATION"])
        wave_radiation = self._parse_Wave_Radiation(self.hydrodynamic["RADIATION_STRESS"])
        sources = self._parse_Sources(self.hydrodynamic["SOURCES"])
        initial_conditions = self._parse_Initial_Condition(self.hydrodynamic["INITIAL_CONDITIONS"])
        temperature_salinity = self._parse_Temperature_Salinity(self.hydrodynamic["TEMPERATURE_SALINITY_MODULE"])
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
                                              "Tidal Potential": tide,
                                              "Precipitation": precipitation,
                                              "Evaporation": evaporation,
                                              "Infiltration": infiltration,
                                              "Wave Radiation": wave_radiation,
                                              "Sources": sources,
                                              "Initial Conditions": initial_conditions,
                                              "Temperature and Salinity": temperature_salinity}
                                            }

    def _parse_Boundaries(self, boundaries):
        self.boundaries = {"CODE_1": {"Name": "Land Boundary"}}
        for key in boundaries.keys():
            if isinstance(boundaries[key], PfsSection):
                self.boundaries = self.boundaries | {key: {"Name": boundaries[key]["name"]}}
        
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

    def _parse_Precipitation(self, precipitation):
        if precipitation["type"] == 0:
            return {"Precipitation": "Not included"}
        elif precipitation["type"] == 1:
            part1 = {"Precipitation type": "Specified precipitation"}
            if precipitation["format"] == 0:
                part1 = part1 | {"Precipitation format": "Constant", "Constant precipitation value": str(precipitation["constant_value"]) + " mm/day"}
            elif precipitation["format"] == 1:
                part1 = part1 | {"Precipitation format": "Varying in time, constant in space"}
            elif precipitation["format"] == 3:
                part1 = part1 | {"Precipitation format": "Varying in time and space"}
            part1 = part1 | {"Soft time interval": str(precipitation["soft_time_interval"]) + " s"}
            return part1
        elif precipitation["type"] == 2:
            part1 = {"Precipitation type": "Net precipitation"}
            if precipitation["format"] == 0:
                part1 = part1 | {"Precipitation format": "Constant", "Constant precipitation value": str(precipitation["constant_value"]) + " mm/day"}
            elif precipitation["format"] == 1:
                part1 = part1 | {"Precipitation format": "Varying in time, constant in space"}
            elif precipitation["format"] == 3:
                part1 = part1 | {"Precipitation format": "Varying in time and space"}
            part1 = part1 | {"Soft time interval": str(precipitation["soft_time_interval"]) + " s"}
            return part1
        
    def _parse_Evaporation(self, evaporation):
        if evaporation["type"] == 0:
            return {"Evaporation": "Not included"}
        elif evaporation["type"] == 1:
            part1 = {"Evaporation type": "Specified evaporation"}
            if evaporation["format"] == 0:
                part1 = part1 | {"Evaporation format": "Constant", "Constant evaporation value": str(evaporation["constant_value"]) + " mm/day"}
            elif evaporation["format"] == 1:
                part1 = part1 | {"Evaporation format": "Varying in time, constant in space"}
            elif evaporation["format"] == 3:
                part1 = part1 | {"Evaporation format": "Varying in time and space"}
            part1 = part1 | {"Soft time interval": str(evaporation["soft_time_interval"]) + " s"}
            return part1

    def _parse_Infiltration(self, infiltration):
        if infiltration["type"] == 0:
            return {"Infiltration": "Not included"}
        elif infiltration["type"] == 1:
            part1 = {"Infiltration type": "Net infiltration rate"}
            if infiltration["format"] == 2:
                part1 = part1 | {"Infiltration format": "Varying in space, constant in time"}
            elif infiltration["format"] == 3:
                part1 = part1 | {"Infiltration format": "Varying in time and space"}
            return part1
        elif infiltration["type"] == 2:
            part1 = {"Infiltration type": "Constant infiltration with capacity", "Infiltration format": "Varying in space, constant in time"}
            return part1
        
    def _parse_Wave_Radiation(self, wave_radiation):
        if wave_radiation["type"] == 0:
            return {"Wave radiation": "Not included"}
        elif wave_radiation["type"] == 1:
            part1 = {"Wave radiation type": "Specified wave radiation", "Wave radiation format": "Varying in time and space", "Soft time interval": str(wave_radiation["soft_time_interval"]) + " s"}
            return part1

    def _parse_Sources(self, sources):
        n_sources = sources["number_of_sources"]
        sources_info = {}
        for i in range(n_sources):
            sources_info[f"Source_{i+1}"] = {}
            source = sources[f"SOURCE_{i+1}"]
            if source["include"] == 0:
                continue
            tmp = {"Source name": source["Name"]}
            location = {"Map projection": self.__parse_coordinates(source["coordinate_type"])}
            if source["interpolation_type"] == 0:
                location = location | {"Coordinates": [source["coordinates"][0], source["coordinates"][1], "Vertical layer number {layer}".format(layer=source["layer"])]}
            elif source["interpolation_type"] == 1:
                location = location | {"Coordinates": [source["coordinates"][0], source["coordinates"][1], "{Z} m below surface".format(Z=source["coordinates"][2])]}
            elif source["interpolation_type"] == 2:
                location = location | {"Coordinates": [source["coordinates"][0], source["coordinates"][1], "{Z} m above bed".format(Z=source["coordinates"][2])]}
            elif source["interpolation_type"] == 3:
                location = location | {"Coordinates": [source["coordinates"][0], source["coordinates"][1], "Z = {Z} m".format(Z=source["coordinates"][2])]}
            tmp = tmp | {"Location": location}
            if source["type"] == 1:
                source_info = {"Source type": "Simple source"}
                if source["format"] == 0:
                    source_info = source_info | {"Source format": "Constant discharge", "Discharge": str(source["constant_value"]) + " m³/s"}
                elif source["format"] == 1:
                    source_info = source_info | {"Source format": "Varying in time"}
                elif source["format"] == 4:
                    source_info = source_info | {"Source format": "Rating curve"}
            elif source["type"] == 2:
                source_info = {"Source type": "Standard source"}
                if source["format"] == 0:
                    source_info = source_info | {"Source format": "Constant discharge",
                                                 "Discharge": str(source["constant_values"][0]) + " m³/s",
                                                 "u-Velocity": str(source["constant_values"][1]) + " m/s",
                                                 "v-Velocity": str(source["constant_values"][2]) + " m/s",
                                                 "w-Velocity": str(source["constant_values"][3]) + " m/s"}
                elif source["format"] == 1:
                    source_info = source_info | {"Source format": "Varying in time"}
            elif source["type"] == 4:
                source_info = {"Source type": "Jet"}
                if source["format"] == 0:
                    source_info = source_info | {"Source format": "Constant discharge", "Discharge": str(source["constant_value"]) + " m³/s"}
                elif source["format"] == 1:
                    source_info = source_info | {"Source format": "Varying in time"}
                source_info = source_info | {"Jet diameter": str(source["diameter"]) + " m",
                                             "Jet horizontal direction angle": str(source["sigma"]) + "°",
                                             "Jet vertical direction angle": str(source["theta"]) + "°",
                                             "Jet maximum travel distance": str(source["maximum_distance"]) + " m"}
                if source["upstream"] == 1:
                    source_info = source_info | {"Jet minimum upstream distance": str(source["distance_upstream"]) + " m"}
            tmp = tmp | {"Source Information": source_info}
            sources_info[f"Source_{i+1}"] = sources_info[f"Source_{i+1}"] | tmp
        if len(sources_info) == 0:
            return {"Sources": "No sources included"}
        else:
            return sources_info
        
    def _parse_Initial_Condition(self, initial_condition):
        if initial_condition["type"] == 1:
            return {"Initial condition type": "Constant initial condition",
                    "Initial surface elevation": str(initial_condition["surface_elevation_constant"]) + " m",
                    "Initial u-velocity": str(initial_condition["u_velocity_constant"]) + " m/s",
                    "Initial v-velocity": str(initial_condition["v_velocity_constant"]) + " m/s",
                    "Initial w-velocity": str(initial_condition["w_velocity_constant"]) + " m/s"}
        elif initial_condition["type"] == 2:
            return {"Initial condition type": "Varying surface elevation in space"}
        elif initial_condition["type"] == 3:
            return {"Initial condition type": "Varying water depth and velocities in space"}

    def _parse_Boundary_Conditions(self, boundary_conditions):
        for key in self.boundaries.keys():
            boundary = boundary_conditions[key]
            if boundary["type"] == 1:
                self.boundaries[key]["Type"] = "Land (zero normal velocity)"
                self.boundaries[key]["Wall friction"] = "Not included" if boundary["type_resistance"] == 0 else str(boundary["resistance_coefficient"]) + " m"
            elif boundary["type"] == 2:
                self.boundaries[key]["Type"] = "Land (zero velocity)"
            elif boundary["type"] == 4:
                self.boundaries[key]["Type"] = "Specified velocity"
                self.__parse_boundary_types(boundary, key, "velocity")
            elif boundary["type"] == 5:
                self.boundaries[key]["Type"] = "Specified flux"
                self.__parse_boundary_types(boundary, key, "flux")
            elif boundary["type"] == 6:
                self.boundaries[key]["Type"] = "Specified water level"
                self.__parse_boundary_types(boundary, key, "level")
            elif boundary["type"] == 7:
                self.boundaries[key]["Type"] = "Specified discharge"
                self.__parse_boundary_types(boundary, key, "discharge")
            elif boundary["type"] == 9:
                self.boundaries[key]["Type"] = "Free outflow"
            elif boundary["type"] == 12:
                self.boundaries[key]["Type"] = "Flather boundary"
                #TODO: Add Flather boundary conditions
                
    def _parse_Temperature_Salinity(self, temperature_salinity):
        part1 = {"Minimum temperature": str(temperature_salinity["EQUATION"]["minimum_temperature"]) + " °C",
                 "Maximum temperature": str(temperature_salinity["EQUATION"]["maximum_temperature"]) + " °C",
                 "Minimum salinity": str(temperature_salinity["EQUATION"]["minimum_salinity"]) + " PSU",
                 "Maximum salinity": str(temperature_salinity["EQUATION"]["maximum_salinity"]) + " PSU"}
        time_integration_scheme = "Low order, fast algorithm" if temperature_salinity["SOLUTION_TECHNIQUE"]["scheme_of_time_integration"] == 1 else "Higher order"
        space_integration_scheme = "Low order, fast algorithm" if temperature_salinity["SOLUTION_TECHNIQUE"]["scheme_of_space_discretization_horizontal"] == 1 else "Higher order"
        part1 = part1 | {"Solution technique": {"Time Integration Scheme": time_integration_scheme, "Space Integration Scheme": space_integration_scheme}}
        part1 = part1 | {"Dispersion": self.__parse_diffusion(temperature_salinity["DIFFUSION"])}
        part1 = part1 | {"Heat Exchange": self.__parse_heat_exchange(temperature_salinity["HEAT_EXCHANGE"])}
        part1 = part1 | {"Precipitation - Evaporation": self.__parse_precipitation_evaporation(temperature_salinity["PRECIPITATION_EVAPORATION"])}
        part1 = part1 | {"Infiltration": self.__parse_infiltration(temperature_salinity["INFILTRATION"])}
        part1 = part1 | {"Sources": self.__parse_sources(temperature_salinity["SOURCES"])}
        return part1

            
                
            

    def __parse_coordinates(self, projection):
        try:
            crs = CRS.from_string(projection)
            coordinate_system = crs.name
        except:
            coordinate_system = projection
        return coordinate_system
    
    def __parse_boundary_types(self, boundary, key, var):
        if var == "velocity":
            vars = ["u-Velocity", "v-Velocity"]
            unit = "m/s"
        elif var == "flux":
            vars = ["p-flux", "q-flux"]
            unit = "m³/s/m"
        if var in ["velocity", "flux"] and boundary["format"] == 0:
            self.boundaries[key]["Format"] = "Constant"
            self.boundaries[key][vars[0]] = str(boundary["constant_values"][0]) + f" {unit}"
            self.boundaries[key][vars[1]] = str(boundary["constant_values"][1]) + f" {unit}"
            self.boundaries[key]["Type of vertical profile"] = "Uniform profile" if boundary["type_of_vertical_profile"] == 1 else "Logarithmic profile"
        elif var in ["level"] and boundary["format"] == 0:
            self.boundaries[key]["Format"] = "Constant"
            self.boundaries[key]["Water Level"] = str(boundary["constant_value"]) + " m"
        elif var in ["velocity", "flux"] and boundary["format"] == 1:
            self.boundaries[key]["Format"] = "Varying in time, constant along boundary"
            self.boundaries[key]["Type of vertical profile"] = "Uniform profile" if boundary["type_of_vertical_profile"] == 1 else "Logarithmic profile"
            self.boundaries[key]["Time interpolation type"] = "Linear" if boundary["type_of_time_interpolation"] == 1 else "Piecewise cubic"
        elif var in ["level"] and boundary["format"] == 1:
            self.boundaries[key]["Format"] = "Varying in time, constant along boundary"
            self.boundaries[key]["Time interpolation type"] = "Linear" if boundary["type_of_time_interpolation"] == 1 else "Piecewise cubic"
        elif var in ["velocity", "flux"] and boundary["format"] == 2:
            self.boundaries[key]["Format"] = "Varying in time and along boundary"
            self.boundaries[key]["Time interpolation type"] = "Linear" if boundary["type_of_time_interpolation"] == 1 else "Piecewise cubic"
            self.boundaries[key]["Space interpolation type"] = "Normal" if boundary["type_of_space_interpolation"] == 1 else "Reverse order"
        elif var in ["level"] and boundary["format"] == 2:
            self.boundaries[key]["Format"] = "Varying in time and along boundary"
            self.boundaries[key]["Time interpolation type"] = "Linear" if boundary["type_of_time_interpolation"] == 1 else "Piecewise cubic"
            self.boundaries[key]["Space interpolation type"] = "Normal" if boundary["type_of_space_interpolation"] == 1 else "Reverse order"
        elif var in ["level"] and boundary["format"] == 4:
            self.boundaries[key]["Format"] = "Rating curve"
        elif var in ["discharge"] and boundary["format"] == 0:
            self.boundaries[key]["Format"] = "Constant"
            self.boundaries[key]["Discharge"] = str(boundary["constant_value"]) + " m³/s"
            self.boundaries[key]["Approach"] = "Weak formulation" if boundary["approach"] == 1 else "Strong formulation"
        elif var in ["discharge"] and boundary["format"] == 1:
            self.boundaries[key]["Format"] = "Varying in time, constant along boundary"
            self.boundaries[key]["Approach"] = "Weak formulation" if boundary["approach"] == 1 else "Strong formulation"
            self.boundaries[key]["Time interpolation type"] = "Linear" if boundary["type_of_time_interpolation"] == 1 else "Piecewise cubic"
            self.boundaries[key]["Type of vertical profile"] = "Uniform profile" if boundary["type_of_vertical_profile"] == 1 else "Logarithmic profile"
        elif var in ["discharge"] and boundary["format"] == 4:
            self.boundaries[key]["Format"] = "Rating curve"
            self.boundaries[key]["Approach"] = "Weak formulation" if boundary["approach"] == 1 else "Strong formulation"
            self.boundaries[key]["Type of vertical profile"] = "Uniform profile" if boundary["type_of_vertical_profile"] == 1 else "Logarithmic profile"
            
        self.boundaries[key]["Type of soft start"] = "Linear variation" if boundary["type_of_soft_start"] == 1 else "Sine variation"
        self.boundaries[key]["Soft time interval"] = str(boundary["soft_time_interval"]) + " s"
        if var in ["velocity", "flux"]:
            self.boundaries[key]["Reference {var}".format(var=vars[0])] = str(boundary["reference_values"][0]) + f" {unit}"
            self.boundaries[key]["Reference {var}".format(var=vars[1])] = str(boundary["reference_values"][1]) + f" {unit}"
        elif var in ["level"]:
            self.boundaries[key]["Reference water level"] = str(boundary["reference_value"]) + " m"
            self.boundaries[key]["Coriolis correction for boundary data"] = "Included" if boundary["type_of_coriolis_correction"] == 1 else "Not included"
            self.boundaries[key]["Wind correction for boundary data"] = "Included" if boundary["type_of_wind_correction"] == 1 else "Not included"
            self.boundaries[key]["Pressure correction for boundary data"] = "Included" if boundary["type_of_pressure_correction"] == 1 else "Not included"
            self.boundaries[key]["Radiation stress correction for boundary data"] = "Included" if boundary["type_of_radiation_stress_correction"] == 1 else "Not included"

    def __parse_diffusion(self, diffusion):
        horizontal = diffusion["HORIZONTAL_DIFFUSION"]
        vertical = diffusion["VERTICAL_DIFFUSION"]
        horizontal_diffusion = self.___parse_dispersion(horizontal)
        vertical_diffusion = self.___parse_dispersion(vertical)
        return {"Horizontal Dispersion": horizontal_diffusion, "Vertical Dispersion": vertical_diffusion}
        
    def __parse_heat_exchange(self, heat_exchange):
        if heat_exchange["type"] == 0:
            output = {"Heat exchange": "Not included"}
        elif heat_exchange["type"] == 1:
            # Latent heat
            latent = {"Conatant in Dalton's law": heat_exchange["Daltons_law_A"],
                   "Wind coefficient in Dalton's law": heat_exchange["Daltons_law_B"],
                   "Critical wind speed": str(heat_exchange["latent_heat_critical_wind_speed"]) + " m/s"}
            # Sensible heat
            sensible = {"Tranfer coefficient for heating": heat_exchange["sensible_heat_transfer_coefficient_heating"],
                        "Transfer coefficient for cooling": heat_exchange["sensible_heat_transfer_coefficient_cooling"],
                        "Critical wind speed": str(heat_exchange["sensible_heat_critical_wind_speed"]) + " m/s"}
            # Short wave radiation
            if heat_exchange["type_of_short_wave_radiation"] == 1:
                short_wave = {"Formulation": "Empirical",
                              "Sun constant, a in Angstroms's law": heat_exchange["Angstroms_law_A"],
                              "Sun constant, b in Angstroms's law": heat_exchange["Angstroms_law_B"],
                              "Displacement (summer time)": str(heat_exchange["displacement_hours"]) + " h",
                              "Standard meridian for time zone": heat_exchange["standard_meridian"],
                              "Beta in Beer's law": heat_exchange["Beers_law_beta"],
                              "Type of abosrption in water column": "Normalized" if heat_exchange["type_of_solar_radiation"] == 1 else "Not normalized"}
            elif heat_exchange["type_of_short_wave_radiation"] == 2:
                short_wave = {"Formulation": "Specified solar radiation",
                              "Displacement (summer time)": str(heat_exchange["displacement_hours"]) + " h",
                              "Standard meridian for time zone": heat_exchange["standard_meridian"],
                              "Beta in Beer's law": heat_exchange["Beers_law_beta"],
                              "Type of abosrption in water column": "Normalized" if heat_exchange["type_of_solar_radiation"] == 1 else "Not normalized"}
                short_wave = short_wave | {"Light extinction coefficient": self.___parse_data(heat_exchange["LIGHT_EXTINCTION"], "light_extinction")}
                short_wave = short_wave | {"Radiation data": self.___parse_data(heat_exchange["SHORT_WAVE_RADIATION_DATA"], "short_wave")}
            elif heat_exchange["type_of_short_wave_radiation"] == 3:
                short_wave = {"Formulation": "Specified net short wave radiation",
                              "Beta in Beer's law": heat_exchange["Beers_law_beta"],
                              "Type of abosrption in water column": "Normalized" if heat_exchange["type_of_solar_radiation"] == 1 else "Not normalized"}
                short_wave = short_wave | {"Light extinction coefficient": self.___parse_data(heat_exchange["LIGHT_EXTINCTION"], "light_extinction")}
                short_wave = short_wave | {"Radiation data": self.___parse_data(heat_exchange["SHORT_WAVE_RADIATION_DATA"], "short_wave")}
            # Long wave radiation
            if heat_exchange["type_of_long_wave_radiation"] == 1:
                long_wave = {"Formulation": "Empirical"}
            elif heat_exchange["type_of_long_wave_radiation"] == 2:
                long_wave = {"Formulation": "Specified atmospheric radiation"}
                long_wave = long_wave | {"Radiation data": self.___parse_data(heat_exchange["LONG_WAVE_RADIATION_DATA"], "long_wave")}
            elif heat_exchange["type_of_long_wave_radiation"] == 3:
                long_wave = {"Formulation": "Specified net long wave radiation"}
                long_wave = long_wave | {"Radiation data": self.___parse_data(heat_exchange["LONG_WAVE_RADIATION_DATA"], "long_wave")}
            # Atmospheric condition
            air_temp = self.___parse_data(heat_exchange["air_temperature"], "air_temp")
            rel_humidity = self.___parse_data(heat_exchange["relative_humidity"], "rel_humidity")
            clear_coef = self.___parse_data(heat_exchange["clearness_coefficient"], "clear_coef")
            if heat_exchange["type_of_short_wave_radiation"] == 1 or heat_exchange["type_of_long_wave_radiation"] == 1:
                atmospheric = air_temp | rel_humidity | clear_coef
            else:
                atmospheric = air_temp | rel_humidity
            # Ground heat
            if heat_exchange["type_of_ground_heat"] == 0:
                ground = "Not included"
            elif heat_exchange["type_of_ground_heat"] == 1:
                ground = {"Distance below ground": str(heat_exchange["distance_below_ground"]) + " m"}
                thermal = self.___parse_data(heat_exchange["THERMAL_CONDUCTIVITY"], "thermal")
                ground_temp = self.___parse_data(heat_exchange["GROUND_TEMPERATURE"], "ground_temp")
                ground = ground | thermal | ground_temp
            output = {"Latent heat": latent, "Sensible heat": sensible, "Short wave radiation": short_wave, "Long wave radiation": long_wave, "Atmospheric condition": atmospheric, "Ground heat": ground}
            
        return output
    
    def __parse_precipitation_evaporation(self, precipitation_evaporation):
        if precipitation_evaporation["type_of_precipitation"] == 1:
            precipitation = {"Precipitation type": "Ambient water temperature"}
        elif precipitation_evaporation["type_of_precipitation"] == 2:
            precipitation = {"Precipitation type": "Specified temperature"} | self.___parse_data(precipitation_evaporation["PRECIPITATION"], "temp", soft_start=True)
        if precipitation_evaporation["type_of_evaporation"] == 1:
            evaporation = {"Evaporation type": "Ambient water temperature"}
        elif precipitation_evaporation["type_of_evaporation"] == 2:
            evaporation = {"Evaporation type": "Specified temperature"} | self.___parse_data(precipitation_evaporation["EVAPORATION"], "temp", soft_start=True)
        return {"Precipitation": precipitation, "Evaporation": evaporation}

    def __parse_infiltration(self, infiltration):
        if infiltration["type_of_infiltration_temperature"] == 1:
            infiltration = {"Infiltration type": "Ambient water temperature"}
        elif infiltration["type_of_infiltration_temperature"] == 2:
            infiltration = {"Infiltration type": "Specified temperature"} | self.___parse_data(infiltration["INFILTRATION"], "temp", soft_start=True)
        return infiltration

    def __parse_sources(self, sources):
        

    def ___parse_dispersion(self, dispersion):
        if dispersion["type"] == 0:
            output = {"Dispersion Type": "No dispersion"}
        elif dispersion["type"] == 1:
            output = {"Dispersion Formulation": "Scaled eddy viscosity formulation"}
            if dispersion["SCALED_EDDY_VISCOSITY"]["format"] == 0:
                output = output | {"Scaled eddy visocsity format": "Constant", "Constant value": dispersion["SCALED_EDDY_VISCOSITY"]["sigma"]}
            elif dispersion["SCALED_EDDY_VISCOSITY"]["format"] == 2:
                output = output | {"Scaled eddy visocsity format": "Varying in space"}
        elif dispersion["type"] == 2:
            output = {"Horizontal Dispersion Formulation": "Dispersion coefficient formulation"}
            if dispersion["DIFFUSION_COEFFICIENT"]["format"] == 0:
                output = output | {"Dispersion coefficient format": "Constant", "Constant value": str(dispersion["DIFFUSION_COEFFICIENT"]["constant_value"]) + " m²/s"}
            elif dispersion["DIFFUSION_COEFFICIENT"]["format"] == 2:
                output = output | {"Dispersion coefficient format": "Varying in space"}
        return output
    
    def ___parse_data(self, radiation, var, soft_start=False):
        units = {"short_wave": "W/m²",
                 "long_wave": "W/m²",
                 "light_extinction": "1/m",
                 "air_temp": "°C",
                 "rel_humidity": "%",
                 "clear_coef": "%",
                 "thermal": "W/m/°C",
                 "ground_temp": "°C",
                 "temp": "°C",}
        full_name = {"short_wave": "Short wave radiation",
                     "long_wave": "Long wave radiation",
                     "light_extinction": "Light extinction",
                     "air_temp": "Air temperature",
                     "rel_humidity": "Relative humidity",
                     "clear_coef": "Clearness coefficient",
                     "thermal": "Thermal conductivity",
                     "ground_temp": "Ground temperature",
                     "temp": "Temperature"}
        if radiation["format"] == 0:
            output = {f"{full_name[var]}": "Constant", f"{full_name[var]} value": str(radiation["constant_value"]) + f" {units[var]}"}
        elif radiation["format"] == 1:
            output = {f"{full_name[var]}": "Varying in time, constant in space"}
        elif radiation["format"] == 3:
            output = {f"{full_name[var]}": "Varying in time and space"}
        if soft_start:
            output = output | {"Soft time interval": str(radiation["soft_time_interval"]) + " s"}
        return output

sample_files = [f for f in os.listdir("sample_input_files")]
fname = os.path.join("sample_input_files", sample_files[9])

pfs = PFS_Parser(fname)
attributes = pfs.parse()
print(attributes)
