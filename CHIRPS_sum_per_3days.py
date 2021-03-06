import arcpy
import os
from datetime import date
from datetime import timedelta
import multiprocessing



#--------------------- environment variable ---------------------#
input_folder = 'C:\\1. DIODINTA\\3. Work\\11. Lain-Lain\\Long Term Flood\\Data\\1982'
output_folder = 'C:\\1. DIODINTA\\3. Work\\11. Lain-Lain\\Long Term Flood\\Data\\3days_calc'
year = '1982'

#--------------------- function definition ----------------------#

def sum_3_days(data_folder, output_folder):
    for i in os.listdir(data_folder):
        if i.endswith(".tif"):
            parseString = i.split('.')
            data_year = parseString[2]
            data_month = parseString[3]
            data_day = parseString[4]
            data_date = date(int(data_year), int(data_month), int(data_day))
            data_2_date = data_date + timedelta(days=1)
            data_3_date = data_date + timedelta(days=2)
            data_file_1 = os.path.join(data_folder, 'chirps-v2.0.'+str(data_date.year)+'.'
                                       +str(data_date.month).zfill(2)+'.'
                                       +str(data_date.day).zfill(2)+'.tif')
            data_file_2 = os.path.join(data_folder, 'chirps-v2.0.' + str(data_2_date.year) + '.'
                                       + str(data_2_date.month).zfill(2) + '.'
                                       + str(data_2_date.day).zfill(2) + '.tif')
            data_file_3 = os.path.join(data_folder, 'chirps-v2.0.' + str(data_3_date.year) + '.'
                                       + str(data_3_date.month).zfill(2) + '.'
                                       + str(data_3_date.day).zfill(2) + '.tif')
            data_file_123_name = 'chirps-v2.0.{0}{1}{2}.3days.tif'.format(str(data_3_date.year),
                                                                            str(data_3_date.month).zfill(2),
                                                                            str(data_3_date.day).zfill(2))
            if os.path.exists(data_file_1) and os.path.exists(data_file_2) and os.path.exists(data_file_3):
                print(str(data_date)+" next 2 days data exist. Start calculating....")
                arcpy.CheckOutExtension("spatial")
                with_null_1 = arcpy.sa.SetNull(data_file_1 < 0, data_file_1)
                with_null_2 = arcpy.sa.SetNull(data_file_2 < 0, data_file_2)
                with_null_3 = arcpy.sa.SetNull(data_file_3 < 0, data_file_3)
                data_123 = arcpy.sa.CellStatistics([with_null_1, with_null_2, with_null_3], "SUM", "DATA")
                data_123.save(os.path.join(output_folder, data_file_123_name))
                print(data_file_123_name+ ' is succesfully created')
                arcpy.CheckInExtension("spatial")
            else:
                print(str(data_date)+" does not have complete 2 following data")

def extract_max(folder_to_extract, output_folder, year):
    listoffile = []
    for data in os.listdir(folder_to_extract):
        if data.endswith(".tif"):
            listoffile.append(os.path.join(folder_to_extract, data))
    print("data to calculate is "+str(len(listoffile)))
    print("start running cell statistics to find maximum rainfall rate in "+year+" .....")
    arcpy.CheckOutExtension("spatial")
    max_data_filename = 'chips-v.2.0.{0}.3days.max.tif'.format(year)
    max_data = arcpy.sa.CellStatistics(listoffile, "MAXIMUM", "DATA")
    max_data.save(os.path.join(output_folder, max_data_filename))
    print(max_data_filename + ' is succesfully created')
    arcpy.CheckInExtension("spatial")


#----------------- calculating chirps 3 days data----------------#
if __name__ == '__main__':
    print("input folder "+input_folder)
    print("output folder "+output_folder)
    folder_3days = os.path.join(output_folder, "calc_3days_"+year)
    os.mkdir(folder_3days)
    print(folder_3days+" is succesfully created")
    print("start calculating.........")
    sum_3_days(input_folder, folder_3days)
    #extract_max(folder_3days, output_folder, year)