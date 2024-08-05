
import os
import pandas as pd
import re
import json

                                            ### Directories and files ###
def load_data_from_directories(main_directory, filename):

    # Create an empty DataFrame to hold the results
    result_df = pd.DataFrame()

    # Loop through each subdirectory in the main directory
    for subdir in os.listdir(main_directory):
        subdir_path = os.path.join(main_directory, subdir)

        # Check if it is a directory
        if os.path.isdir(subdir_path):
            file_path = os.path.join(subdir_path, filename)

            # Check if the file exists in the subdirectory
            if os.path.isfile(file_path):
                # Read the file into a DataFrame and append to the list
                df = pd.read_csv(file_path)
                # Get the last row of the DataFrame
                last_row = df.iloc[[-1]]

                # Rename the row index to the filename (without extension)
                row_name = os.path.splitext(subdir)[0]
                last_row.index = [row_name]

                # Append the last row to the result DataFrame
                result_df = pd.concat([result_df, last_row])


    return result_df

def process_files(main_directory, filename):
    # Create an empty DataFrame to hold the results
    combined_df = pd.DataFrame()
    columns_to_analyse = ['Component', 'Component_type', 'Component_code', 'Capacity_kW']
    # Loop through each subdirectory in the main directory
    for subdir in os.listdir(main_directory):
        if subdir == 'current_DES' or subdir == 'debugging':
            continue
        subdir_path = os.path.join(main_directory, subdir)

        # Check if it is a directory
        if os.path.isdir(subdir_path):
            file_path = os.path.join(subdir_path, filename)

            # Check if the file exists in the subdirectory
            if os.path.isfile(file_path):
                # Read the file into a DataFrame
                df = pd.read_csv(file_path)
                df = df[columns_to_analyse]

                # Extract the numerical part of the folder name
                match = re.search(r'\d+', subdir)
                if match:
                    num_part = match.group(0)

                    # Rename the columns by appending the numerical part
                    df.columns = [f"{col}_{num_part}" if col != 'Supply_System' else col for col in df.columns]

                # Append the DataFrame to the combined DataFrame
                combined_df = pd.concat([combined_df, df], axis=1)

    return combined_df

def calculate_percentage_change(df):

    supply_systems = df['Supply_System']
    base_row = df.loc['current_DES']
    percentage_df = (df.iloc[:, 1:] - base_row[1:]) / base_row[1:] * 100
    percentage_df = percentage_df.apply(pd.to_numeric, errors='coerce')
    percentage_df.insert(0, 'Supply_System', supply_systems)

    return percentage_df

def combine_selected_system_with_structure(selected_systems, structure_df):
    # Combine the selected systems with the structure DataFrame
    selected_systems_df = pd.DataFrame(selected_systems, columns=structure_df.columns)
    combined_df = pd.concat([selected_systems_df, structure_df], axis=1)

    return combined_df

def process_energy_system_data(main_directory, selected_systems, filename_structure, selected_systems_structure, context, scenario, dict_availabilities):
    """
    Process energy system data from CSV files and update the growing DataFrame.

    Parameters:
    - base_path: A string representing the fixed part of the path.
    - variable_parts: A list of strings representing the variable parts of the path.
    - growing_df: A DataFrame to which the extracted data will be appended.

    Returns:
    - The updated DataFrame.
    """
    system_names = selected_systems['Supply_System']

    # Construct the path
    for name in system_names:
        file_path = os.path.join(main_directory, name, filename_structure)

        if not os.path.isfile(file_path):
            continue

        # Read the CSV file
        df = pd.read_csv(file_path)

        # Extract relevant columns
        components = df['Component']
        capacities = df['Capacity_kW'].values
        codes = df['Component_code']
        combined_components_code = [f"{code}_{component}" for code, component in zip(codes, components)]

        temp_df = pd.DataFrame(columns=combined_components_code)
        temp_df.loc[0] = capacities

        # Drop duplicate columns, keep first occurrence
        temp_df = temp_df.loc[:, ~temp_df.columns.duplicated(keep='first')]

        if 'System_name' not in temp_df.columns:
            temp_df.insert(0, 'System_name', name)
            temp_df.insert(1, 'Scenario', scenario)
            temp_df.insert(2, 'Availability', dict_availabilities[context])

        else:
            temp_df['System_name'] = name
            temp_df['Scenario'] = scenario
            temp_df['Availability'] = dict_availabilities[context]

        # Create a DataFrame with the components as columns and the system name as the index
        temp_df = merge_pv_columns(temp_df)
        if selected_systems_structure.empty:
            selected_systems_structure = temp_df
        else:
            selected_systems_structure = pd.concat([selected_systems_structure, temp_df], axis=0, ignore_index=True)

    return selected_systems_structure

def line_plot_dataframe_preprocessing(combined_df, hours_to_plot):

    # Clean up not useful columns
    if 'PV1_E230AC_output' in combined_df.columns:
        combined_df = combined_df.drop(columns='PV1_E230AC_output')
    if 'PV2_E230AC_output' in combined_df.columns:
        combined_df = combined_df.drop(columns='PV2_E230AC_output')
    if 'PV3_E230AC_output' in combined_df.columns:
        combined_df = combined_df.drop(columns='PV3_E230AC_output')

    # Select only 24 hours of operation and exclude insignificant values
    carriers_to_plot = combined_df.iloc[:hours_to_plot, :]
    columns_to_drop = carriers_to_plot.columns[(carriers_to_plot < 0.5).all()]
    carriers_to_plot = carriers_to_plot.drop(columns=columns_to_drop)

    # Unify energy demand from auxiliary systems
    for column in carriers_to_plot.columns:
        if 'E230AC' in column:
            if 'input' in column:
                if 'E230AC_demand' in carriers_to_plot.columns:
                    carriers_to_plot['E230AC_demand'] += carriers_to_plot[column]
                else:
                    carriers_to_plot['E230AC_demand'] = carriers_to_plot[column]
                carriers_to_plot = carriers_to_plot.drop(columns=column)

        if 'main' in column:
            parts = column.split('_')
            if 'SC' in parts[0] or 'PV' in parts[0] or 'CH' in parts[0] or 'FU' in parts[0]:
                carriers_to_plot = carriers_to_plot.rename(columns={column: f'{parts[0]}_{parts[1]}_output'})
            elif 'BT' in parts[0] or 'TES' in parts[0]:
                carriers_to_plot = carriers_to_plot.rename(columns={column: f'{parts[0]}_{parts[1]}_SOC'})

        if 'T30W_main' in column:
            carriers_to_plot = carriers_to_plot.drop(columns=column)

    return carriers_to_plot

def combine_excel_sheets(file_path):
    excel_data = pd.read_excel(file_path, sheet_name=None)
    combined_df = pd.DataFrame()

    for sheet_name, sheet_data in excel_data.items():
        sheet_data = sheet_data.iloc[:, 1:]  # Drop the indexing column
        if sheet_name == 'bought_carriers' or sheet_name == 'sold_carriers':
            parts = sheet_name.split('_')
            sheet_data = sheet_data.rename(columns={col: f'{col}_{parts[0]}' for col in sheet_data.columns})
        combined_df = pd.concat([combined_df, sheet_data], axis=1)

    return combined_df

def compute_mean_and_std(df):
    mean_df = pd.DataFrame()
    std_df = pd.DataFrame()

    for column in df.columns:
        profile = df[column]
        daily_profiles = []

        for i in range(365):
            daily_profile = profile[i*24:(i+1)*24].reset_index(drop=True)
            daily_profiles.append(daily_profile)

        daily_profiles_df = pd.concat(daily_profiles, axis=1)
        mean_df[column] = daily_profiles_df.mean(axis=1)
        std_df[column] = daily_profiles_df.std(ddof=0, axis=1)

    return mean_df, std_df

def merge_pv_columns(df):
    # Check if PV1, PV2, PV3 columns exist
    pv_columns = ['PV1_Photovoltaic Panels', 'PV2_Photovoltaic Panels', 'PV3_Photovoltaic Panels']
    existing_pv_columns = [col for col in pv_columns if col in df.columns]

    # Create a new column 'PV' with summed capacities
    df['PV_Photovoltaic Panels'] = df[existing_pv_columns].sum(axis=1)

    # Drop the old PV columns
    df = df.drop(columns=existing_pv_columns)

    return df
def get_connected_buildings(geojson_path):
    # Load the GeoJSON file
    with open(geojson_path, 'r') as file:
        data = json.load(file)

    # Extract building names from the GeoJSON file
    buildings = []
    for feature in data['features']:
        if 'properties' in feature and 'Building' in feature['properties']:
            building_id = feature['properties']['Building']
            if building_id != None and building_id.startswith('B10'):
                buildings.append(building_id)
    buildings = sorted(set(buildings))

    return buildings


def get_district_buildings(demand_folder_path):
    # List all files in the demand folder
    demand_files = os.listdir(demand_folder_path)

    # Extract building names from the filenames
    buildings = []
    for filename in demand_files:
        if filename.startswith('B'):
            building_name = filename.split('.')[0]
            buildings.append(building_name)

    # Get unique and sorted building names
    buildings = sorted(set(buildings))

    return buildings

def process_scenario(connectivity_df, folder, directory_to_file, geojson, demand, df_scenarios, scenario_name):

    # Filter the dataframe for the given scenario
    filtered_df = df_scenarios[df_scenarios['Scenario'] == scenario_name]

    # Create a list of system names from the filtered dataframe
    system_names = filtered_df['Supply_System'].tolist()

    # Create a list of components availability from the filtered dataframe
    availabilities = filtered_df['Availability'].tolist()

    # Initialize a list to store results
    results = []

    # Construct the path to the demand folder
    demand_folder_path = os.path.join(folder, 'THESIS_TEST_CASES_BASE', scenario_name, demand)

    # Get the list of buildings in the district
    district_buildings = get_district_buildings(demand_folder_path)

    # Iterate over each system
    for j, system in enumerate(system_names):
        # Get the availability context
        if availabilities[j] == 'No_renewables':
            context = 'THESIS_TEST_CASES_BASE'
        else:
            context = 'THESIS_TEST_CASES_RENEWABLES'

        # Construct the path to the geojson file
        base_path = os.path.join(folder, context, scenario_name)
        geojson_path =  base_path + directory_to_file + '/' + system + '/' + geojson

        # Get the list of connected buildings
        connected_buildings = get_connected_buildings(geojson_path)

        # Calculate the percentage of connected buildings
        connected_count = len(set(connected_buildings) & set(district_buildings))
        total_buildings = len(district_buildings)
        percentage_connected = (connected_count / total_buildings) * 100 if total_buildings > 0 else 0

        # Store the result in a dataframe
        results.append({
            'scenario': scenario_name,
            'system': system,
            'availability': availabilities[j],
            'connected_percentage': percentage_connected,
            'connected_buildings': connected_buildings,
            'district_buildings': district_buildings
        })

    # Convert results to a dataframe
    if connectivity_df.empty:
        connectivity_df = pd.DataFrame(results)
    else:
        connectivity_df = pd.concat([connectivity_df, pd.DataFrame(results)])

    return connectivity_df
