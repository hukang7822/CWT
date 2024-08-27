'''
Python code to calculate PSCF and CWT analysis for time series
From: Kang Hu, Nanjing University of Information Science and Technology, Hong Liao group
Author: Kang Hu, NUIST
Written: 2024-04-26
Version 1.0
email: 200060@nuist.edu.cn
First column of PM data is date; Second column is PM values
Tdump file and PM file should be saved in different folders
'''


import pysplit
import pandas as pd
import numpy as np
import sys
import os
import math
from datetime import datetime
from pytz import timezone

def File_load_HYSPLT_append(file_path, date_and_time_check, temp):
    data = open(file_path, 'r',encoding='gb18030' , errors='ignore')
    val = data.readline()
    val_split = val.strip().split(' ')
    for i in val_split:
        if len(i) == 0:
            continue
        else:
            val = i
            break
    for i in range(int(val) + 29):
        data.readline()
    while 1:
        line = data.readline()
        if not line:
            break
        line_split = line.strip().split(' ')
        for i in line_split:
            if len(i) == 0:
                continue
            else:
                temp.append(i)
    temp = pd.Series(temp)
    raw_num = len(temp) // 13
    matrix = temp.values.reshape(raw_num, 13)
    Lat_temp = matrix[:, 9]
    Lon_temp = matrix[:, 10]
    kk = 0
    latitude = []
    longitude = []
    while kk < 27:
        try:
            for i in range(25):
                latitude.append(Lat_temp[i*27+kk])
                longitude.append(Lon_temp[i*27+kk])
        except IndexError:
            print(date_and_time_check)
        kk += 1
    latitude = pd.Series(latitude)
    longitude = pd.Series(longitude)
    lat_num = len(latitude) // 25
    matrix_lat = latitude.values.reshape(lat_num, 25)
    matrix_lon = longitude.values.reshape(lat_num, 25)
    return matrix_lat,matrix_lon


def Pollution_read(pwd_path_PM):
    pwd = pwd_path_PM
    data = open(pwd, 'r')
    data.readline()
    time_PM = []
    PM = []
    while 1:
        line = data.readline()
        if not line:
            break
        line_split = line.strip().split(',')
        Date_Time = line_split[0]
        Year = Date_Time.split(' ')[0].split('/')[0]
        Month =  Date_Time.split(' ')[0].split('/')[1]
        Day = Date_Time.split(' ')[0].split('/')[2]
        if int(Month) < 10 and len(Month) < 2:
            Month = '0' + str(Month)
        if int(Day) < 10 and len(Day) < 2:
            Day = '0' + str(Day)
        time_res = Year + '/' + Month + '/' + Day + ' ' + Date_Time.split(' ')[1]
        time_PM.append(time_res)
        PM.append(line_split[1])
    time_PM = pd.Series(time_PM).astype(str)
    PM = pd.Series(PM).astype(float)
    return time_PM, PM

def SameTime_Pick(time_PM,PM,time_HYSPLIT,latitude, longitude, nn, output_file_path):
    pwd = output_file_path
    num2 = len(time_PM)
    while nn < num2:
        time1 = datetime.strptime(time_PM[nn], "%Y/%m/%d %H:%M:%S")
        time2 = datetime.strptime(time_HYSPLIT, "%Y/%m/%d %H:%M:%S")
        if time1 < time2:
            nn += 1
            if nn >= num2:
                break
        elif time1 > time2:
            nn -= 1
            break
        else:
            for kk in range(27):
                res_lat = latitude[kk,:]
                res_lon = longitude[kk,:]
                res_lat = pd.Series(res_lat).astype(float)
                res_lat = pd.concat([res_lat], axis=1)
                res_lon = pd.Series(res_lon).astype(float)
                res_lon = pd.concat([res_lon], axis=1)
                if kk == 0:
                    result = pd.concat([res_lat, res_lon], axis=1)
                else:
                    result = pd.concat([result, res_lat, res_lon], axis=1)
            name_temp = str(time_PM[nn]).split(' ')[0].split('/')
            output_path = pwd + '/' +  'MatrixofTrajPM-' + name_temp[0]
            if os.path.exists(output_path):
                print(name_temp[0] + '-' + name_temp[1] + '-' + name_temp[2])
            else:
                os.mkdir(output_path)
                print(name_temp[0] + '-' + name_temp[1] + '-' + name_temp[2])
            name = output_path + '/' + name_temp[0] + name_temp[1] + name_temp[2] + str(time_PM[nn]).split(' ')[1].split(':')[0] + '-' + str(PM[nn]) + '-PM' + '.txt'
            result.to_csv(name, sep='\t', index=False, header=False)
            break
    return nn


def Grid_points_passed_by_two_coordinates(matrix, FirstLat, FirstLon, SecondLat, SecondLon,PM_value):

    Left_lon_temp, Right_lon_temp, Up_lat_temp, Down_lat_temp, step, weight = Grid_set_up()

    if FirstLat < SecondLat and (FirstLat + step) > SecondLat :
        SecondLat = FirstLat
    elif FirstLat > SecondLat and (FirstLat - step) < SecondLat:
        SecondLat = FirstLat
    if FirstLon < SecondLon and (FirstLon + step) > SecondLon :
        SecondLon = FirstLon
    elif FirstLon > SecondLon and (FirstLon - step) < SecondLon:
        SecondLon = FirstLon
    if FirstLat < Down_lat_temp or FirstLat > Up_lat_temp or SecondLat < Down_lat_temp or SecondLat > Up_lat_temp:
        print ('Error: HYSPLIT trajectory is out of Grid boundary')
        sys.exit()
    if FirstLon < Left_lon_temp or FirstLon > Right_lon_temp or SecondLon < Left_lon_temp or SecondLon > Right_lon_temp:
        print ('Error: HYSPLIT trajectory is out of Grid boundary')
        sys.exit()

    '''
    转换成Grid里的点对应的经纬度
    '''
    startLon_grid,endLon_grid,startLat_grid,endLat_grid = LocationDecision(FirstLat, FirstLon, SecondLat, SecondLon)      ####  经纬度范围 + 1附近点位

    grid_points = set()

    '''
    if startLon_grid > endLon_grid:
        tempVal = endLon_grid
        endLon_grid = startLon_grid
        startLon_grid = tempVal
    if startLat_grid > endLat_grid:
        tempVal = endLat_grid
        endLat_grid = startLat_grid
        startLat_grid = tempVal
    '''

    if endLon_grid == startLon_grid:
        for input_grid in range(int(startLat_grid), int(endLat_grid) + 1, 1):
            grid_points = line_through_grid_y_append(input_grid, startLon_grid, grid_points)
    elif startLon_grid != endLon_grid and startLat_grid != endLat_grid:
        slope = (float(SecondLat) - float(FirstLat)) / (float(SecondLon) - float(FirstLon))
        intercep = FirstLat - FirstLon * slope
        for input_grid in range(startLon_grid, endLon_grid + 1):
            grid_points = line_through_grid_x(input_grid, step, Left_lon_temp, Down_lat_temp, slope, intercep,
                                              startLat_grid, endLat_grid, grid_points)
        for input_grid in range(startLat_grid, endLat_grid + 1):
            grid_points = line_through_grid_y(input_grid, step, Left_lon_temp, Down_lat_temp, slope, intercep,
                                              startLon_grid, endLon_grid, grid_points)
    if startLat_grid == endLat_grid :
        for input_grid in range(startLon_grid, endLon_grid + 1):
            grid_points = line_through_grid_x_append(input_grid, startLat_grid, grid_points)

    elif startLat_grid != endLat_grid and startLon_grid != endLon_grid:
        slope = (float(SecondLat) - float(FirstLat)) / (float(SecondLon) - float(FirstLon))
        intercep = FirstLat - FirstLon * slope
        for input_grid in range(startLon_grid, endLon_grid + 1):
            grid_points = line_through_grid_x(input_grid, step, Left_lon_temp, Down_lat_temp, slope, intercep, startLat_grid, endLat_grid, grid_points)
        for input_grid in range(startLat_grid, endLat_grid + 1):
            grid_points = line_through_grid_y(input_grid, step, Left_lon_temp, Down_lat_temp, slope, intercep, startLon_grid, endLon_grid, grid_points)
    for values in grid_points:
        matrix = Grid_val(matrix, values[1], values[0], PM_value)
    return matrix


    """
    FirstLat_point = (float(FirstLat) - float(Down_lat_temp)) / float(step)
    FirstLon_point = (float(FirstLon) - float(Left_lon_temp)) / float(step)
    SecondLat_point = (float(SecondLat) - float(Down_lat_temp)) / float(step)
    SecondLon_point = (float(SecondLon) - float(Left_lon_temp)) / float(step)
    start_point = [FirstLat_point, SecondLat_point]
    second_point = [FirstLon_point, SecondLon_point]
    #plt.xticks(np.arange(0,15,1))
    #plt.yticks(np.arange(0,15,1))
    plt.plot(start_point,second_point, color = 'b')
    plt.imshow(matrix.T, cmap = plt.cm.hot, vmin = 0, vmax = 10)
    #plt.contourf(arra_lon, arra_lat, matrix, levels = 10, colors = 'r')
    plt.grid(True)
    """

def FindNearstLessThan(down, up, step, inputVar):
    num = int((float(up) - float(down)) / float(step)) + 1
    for i in range(num):
        asdfasf = float(down) + float(step) * int(i)
        if float(down) + float(step) * int(i) == float(inputVar):
            return i
        elif float(down) + float(step) * int(i) < float(inputVar) and float(down) + float(step) * int(i+1) > float(inputVar):
            return i

def FindNearstMoreThan(down, up, step, inputVar):
    num = int((float(up) - float(down)) / float(step)) + 1
    for i in range(num):
        if float(down) + float(step) * int(i) == float(inputVar):
            return i
        elif float(down) + float(step) * int(i) < float(inputVar) and float(down) + float(step) * int(i + 1) > float(
                inputVar):
            return i+1

def LocationDecision(FirstLat, FirstLon, SecondLat, SecondLon):
    Left_lon_temp, Right_lon_temp, Up_lat_temp, Down_lat_temp, step, weight = Grid_set_up()

    if FirstLat == SecondLat:
        leftLat_boundary = FindNearstLessThan(Down_lat_temp, Up_lat_temp, step, FirstLat)
        rightLat_boundary = leftLat_boundary
        if FirstLon < SecondLon:
            leftLon_boundary = FindNearstMoreThan(Left_lon_temp, Right_lon_temp, step, FirstLon)
            rightLon_boundary = FindNearstLessThan(Left_lon_temp, Right_lon_temp, step, SecondLon)
        elif FirstLon > SecondLon:
            leftLon_boundary = FindNearstMoreThan(Left_lon_temp, Right_lon_temp, step, SecondLon)
            rightLon_boundary = FindNearstLessThan(Left_lon_temp, Right_lon_temp, step, FirstLon)
    if FirstLon == SecondLon:
        leftLon_boundary = FindNearstLessThan(Left_lon_temp, Right_lon_temp, step, FirstLon)
        rightLon_boundary = leftLon_boundary
        if FirstLat < SecondLat:
            leftLat_boundary = FindNearstMoreThan(Down_lat_temp, Up_lat_temp, step, FirstLat)
            rightLat_boundary = FindNearstLessThan(Down_lat_temp, Up_lat_temp, step, SecondLat)
        elif FirstLat > SecondLat:
            leftLat_boundary = FindNearstMoreThan(Down_lat_temp, Up_lat_temp, step, SecondLat)
            rightLat_boundary = FindNearstLessThan(Down_lat_temp, Up_lat_temp, step, FirstLat)

    if FirstLon < SecondLon and FirstLat < SecondLat:
        leftLon_boundary = FindNearstMoreThan(Left_lon_temp, Right_lon_temp, step, FirstLon)
        rightLon_boundary = FindNearstLessThan(Left_lon_temp, Right_lon_temp, step, SecondLon)
        leftLat_boundary = FindNearstMoreThan(Down_lat_temp, Up_lat_temp, step, FirstLat)
        rightLat_boundary = FindNearstLessThan(Down_lat_temp, Up_lat_temp, step, SecondLat)
    elif FirstLon < SecondLon and FirstLat > SecondLat:
        leftLon_boundary = FindNearstMoreThan(Left_lon_temp, Right_lon_temp, step, FirstLon)
        rightLon_boundary = FindNearstLessThan(Left_lon_temp, Right_lon_temp, step, SecondLon)
        leftLat_boundary = FindNearstMoreThan(Down_lat_temp, Up_lat_temp, step, SecondLat)
        rightLat_boundary = FindNearstLessThan(Down_lat_temp, Up_lat_temp, step, FirstLat)
    elif FirstLon > SecondLon and FirstLat < SecondLat:
        leftLon_boundary = FindNearstMoreThan(Left_lon_temp, Right_lon_temp, step, SecondLon)
        rightLon_boundary = FindNearstLessThan(Left_lon_temp, Right_lon_temp, step, FirstLon)
        leftLat_boundary = FindNearstMoreThan(Down_lat_temp, Up_lat_temp, step, FirstLat)
        rightLat_boundary = FindNearstLessThan(Down_lat_temp, Up_lat_temp, step, SecondLat)
    elif FirstLon > SecondLon and FirstLat > SecondLat:
        leftLon_boundary = FindNearstMoreThan(Left_lon_temp, Right_lon_temp, step, SecondLon)
        rightLon_boundary = FindNearstLessThan(Left_lon_temp, Right_lon_temp, step, FirstLon)
        leftLat_boundary = FindNearstMoreThan(Down_lat_temp, Up_lat_temp, step, SecondLat)
        rightLat_boundary = FindNearstLessThan(Down_lat_temp, Up_lat_temp, step, FirstLat)
    return leftLon_boundary,rightLon_boundary,leftLat_boundary,rightLat_boundary

def Point2Value(point, step, startVar):
    res = float(point) * float(step) + float(startVar)
    return res

def line_through_grid_x(input_grid, step, Left_lon_temp, Down_lat_temp, slope, intercep, startLat_grid, endLat_grid, grid_points):
    input_val = Point2Value(input_grid, step, Left_lon_temp)
    y_val = float(input_val) * float(slope) + float(intercep)
    y_val = round(y_val, 2)
    for i in range(startLat_grid, endLat_grid + 1, 1):
        temp_val = Point2Value(i, step, Down_lat_temp)
        temp_val_next = Point2Value(i + 1, step, Down_lat_temp)
        if float(y_val) == float(temp_val):
            grid_points.add((input_grid - 1, i - 1))
            grid_points.add((input_grid - 1, i))
            grid_points.add((input_grid, i - 1))
            grid_points.add((input_grid, i))
            break
        elif float(y_val) > float(temp_val) and float(y_val) < float(temp_val_next):
            grid_points.add((input_grid - 1, i))
            grid_points.add((input_grid, i))
            break
    return grid_points

def line_through_grid_y(input_grid, step, Left_lon_temp, Down_lat_temp, slope, intercep, startLon_grid, endLon_grid, grid_points):
    input_val = Point2Value(input_grid, step, Down_lat_temp)
    x_val = (float(input_val) - float(intercep)) / float(slope)
    x_val = round(x_val,2)
    for i in range(startLon_grid, endLon_grid + 1, 1):
        temp_val = Point2Value(i, step, Left_lon_temp)
        temp_val_next = Point2Value(i + 1, step, Left_lon_temp)
        if float(x_val) == float(temp_val):
            grid_points.add((i - 1, input_grid - 1))
            grid_points.add((i, input_grid - 1))
            grid_points.add((i - 1, input_grid))
            grid_points.add((i, input_grid))
            break
        elif float(x_val) > float(temp_val) and float(x_val) < float(temp_val_next):
            grid_points.add((i, input_grid - 1))
            grid_points.add((i, input_grid))
            break
    return grid_points

def line_through_grid_y_append(input_grid, startLon_grid, grid_points):
    grid_points.add((startLon_grid - 1, input_grid + 1))
    grid_points.add((startLon_grid - 1, input_grid))
    grid_points.add((startLon_grid, input_grid + 1))
    grid_points.add((startLon_grid, input_grid))
    return grid_points

def line_through_grid_x_append(input_grid, startLat_grid, grid_points):
    grid_points.add((input_grid - 1, startLat_grid - 1))
    grid_points.add((input_grid - 1, startLat_grid))
    grid_points.add((input_grid, startLat_grid - 1))
    grid_points.add((input_grid, startLat_grid))
    return grid_points

def Grid_val(matrix, value1, value2, PM_value):
    value1 = int(value1)
    value2 = int(value2)
    if math.isnan(matrix[value1][value2]):
        matrix[value1][value2] = float(PM_value)
    else:
        matrix[value1][value2] += float(PM_value)
    return matrix

def CWT_cal(file_path, numerator, denominator):
    data = pd.read_csv(file_path, sep='\t', encoding='ISO-8859-1', error_bad_lines=False)
    pm_var = float(str(file_path).split('-')[2])
    for j in range(27):
        lat_var = data.iloc[:, 2*j]
        lon_var = data.iloc[:, 2*j+1]
        for i in range(data.shape[0] - 1):
        #for i in range(26, 29, 1):
            numerator = Grid_points_passed_by_two_coordinates(numerator, lat_var[i], lon_var[i], lat_var[i + 1], lon_var[i + 1], pm_var)
            denominator = Grid_points_passed_by_two_coordinates(denominator, lat_var[i], lon_var[i], lat_var[i + 1], lon_var[i + 1], 1)

    return numerator, denominator

def CWT_append(numerator, denominator, W, Weight):
    denominator_mean = np.nanmean(denominator)
    condition1 = (denominator > 3 * denominator_mean)
    W[condition1] = Weight[0]
    condition2 = (denominator > 1.5 * denominator_mean) & (denominator <= 3 * denominator_mean)
    W[condition2] = Weight[1]
    condition3 = (denominator > denominator_mean) & (denominator <= 1.5 * denominator_mean)
    W[condition3] = Weight[2]
    condition4 = (denominator <= denominator_mean)
    W[condition4] = Weight[3]
    res = numerator / denominator * W
    return res

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def HYSPLTI():
    working_dir = '/Volumes/Macintosh HD - Data/Users/zju/hysplit/working'
    pwd = '/Volumes/HK/NUIST/CWT/HYSPLIT/'
    meteo_dir = '/Volumes/Hukang/gdas/'
    basename = 'traj-test'
    for i in range(2013,2021):
        storage_dir = pwd + str(i)
        years = [i]
        months = [1,2,3,4,5,6,7,8,9,10,11,12]
        hours = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]
        altitudes = [50]
        location = (39.98, 116.33)
        runtime = -24

        pysplit.generate_bulktraj(basename, working_dir, storage_dir, meteo_dir,
                                  years, months, hours, altitudes, location, runtime,
                                  monthslice=slice(0, 32, 1), get_reverse=False,
                                  get_clipped=False, hysplit='/Users/zju/hysplit/exec/hyts_ens')

def UTC_to_BJ_Time():
    pwd = '/Volumes/HK/NUIST/CWT/HYSPLIT/2020'

    for root, dirs, files in os.walk(pwd):
        files[:] = [d for d in files if not d.startswith('.')]
        files.sort(reverse=True)
        for file in files:
            file_path = os.path.join(root, file)
            UTC_temp = file_path.split('/')[-1][-10:]
            name_temp = file_path.split('/')[-1][0:9]
            year = UTC_temp[:4]
            month = UTC_temp[4:6]
            day = UTC_temp[6:8]
            hour = UTC_temp[8:10]
            UTC = datetime(year=int(year),month = int(month),day = int(day),hour = int(hour), tzinfo=timezone('UTC'))
            BJ_time_temp = UTC.astimezone(timezone('Asia/Shanghai'))
            print (BJ_time_temp)
            BJ_time_temp1 = str(BJ_time_temp).split('+')[0]
            BJ_time_temp2 = str(BJ_time_temp1).split(' ')[0].split('-')
            BJ_time = BJ_time_temp2[0] + BJ_time_temp2[1] + BJ_time_temp2[2] + str(BJ_time_temp1).split(' ')[1][0:2]
            name_new = pwd + '/' + name_temp + '-' + BJ_time
            os.rename(file_path, name_new)

#####输出时间相同的HYSPLIT和PM数据在MatrixofTrajPM-Year文件夹中
def File_load_HYSPLT():
    import os
    pwd_tdump = '/Volumes/HK/NUIST/CWT/HYSPLIT/2019/'
    pwd_PM = '/Volumes/HK/NUIST/CWT/PM_data/PM25.csv'
    pwd_output = '/Volumes/HK/NUIST/CWT'
    nn = 0
    time_PM, PM = Pollution_read(pwd_PM)
    for root, dirs, files in os.walk(pwd_tdump):
        files[:] = [d for d in files if not d.startswith('.')]
        files.sort()
        for file in files:
            year_ = file[-10:-6]
            month_ = file[-6:-4]
            day_ = file[-4:-2]
            hour_ = file[-2:]
            date_and_time_check = year_ + '/' + month_ + '/' + day_ + ' ' + hour_ + ':00:00'
            file_path = os.path.join(root, file)
            temp = []
            matrix_lat,matrix_lon = File_load_HYSPLT_append(file_path, date_and_time_check, temp)
            nn = SameTime_Pick(time_PM, PM, date_and_time_check, matrix_lat,matrix_lon, nn, pwd_output)
            nn += 1

def Grid_set_up():
    Left_lon = -180
    Right_lon = 180
    Up_lat = 90
    Down_lat = -90
    Degree = 0.25
    Weight = [1, 0.7, 0.4, 0.17]
    return Left_lon, Right_lon, Up_lat, Down_lat, Degree, Weight

def CWT():

    pwd = '/Volumes/HK/NUIST/CWT/append-a'
    f = '/Volumes/HK/NUIST/CWT/'

    Left_lon_temp, Right_lon_temp, Up_lat_temp, Down_lat_temp, step, Weight = Grid_set_up()
    num_lat = abs(Up_lat_temp - Down_lat_temp) / step + 1
    num_lon = abs(Right_lon_temp - Left_lon_temp) / step + 1
    numerator = np.full((math.ceil(num_lat), math.ceil(num_lon)), np.nan)
    W = np.full((math.ceil(num_lat), math.ceil(num_lon)), np.nan)
    denominator = np.full((math.ceil(num_lat), math.ceil(num_lon)), np.nan)

    for root, dirs, files in os.walk(pwd):
        files[:] = [d for d in files if not d.startswith('.')]
        files.sort()
        for file in files:
            file_path = os.path.join(root, file)
            date_and_time = file[:10]

            numerator, denominator = CWT_cal(file_path, numerator, denominator)
            res = CWT_append(numerator, denominator, W, Weight)

            output_path = f + 'CWT/' + file[:4]
            output = f + 'CWT/' + file[:4] + '/CWT-' + date_and_time + '.csv'
            if os.path.exists(output_path) == False:
                os.mkdir(output_path)
            res = pd.DataFrame(res)
            res.to_csv(output, sep=',')
            print(file)

            numerator = np.full((math.ceil(num_lat), math.ceil(num_lon)), np.nan)
            W = np.full((math.ceil(num_lat), math.ceil(num_lon)), np.nan)
            denominator = np.full((math.ceil(num_lat), math.ceil(num_lon)), np.nan)


def Region_seperate():
    Left_lon_temp, Right_lon_temp, Up_lat_temp, Down_lat_temp, step, Weight = Grid_set_up()
    start_lon_val = 108
    end_lon_val = 121
    start_lat_val = 34
    end_lat_val = 43
    start_lat = math.ceil((start_lat_val+90)/step)
    start_lon = math.ceil((start_lon_val+180)/step)
    end_lat = math.ceil((end_lat_val+90)/step)
    end_lon = math.ceil((end_lon_val+180)/step)
    result = []
    pwd = '/Volumes/HK/NUIST/CWT/CWT/2013'

    for root, dirs, files in os.walk(pwd):
        files[:] = [d for d in files if not d.startswith('.')]
        files.sort()
        for file in files:
            file_path = os.path.join(root, file)
            North = 0
            West = 0
            South = 0
            East = 0
            local = 0
            data = pd.read_csv(file_path, sep=',', encoding='ISO-8859-1', error_bad_lines=False, header = None)
            date_ = file[4:14]      ####
            date = date_[:4] +  '/' + date_[4:6] + '/' + date_[6:8] + ' ' + date_[8:] + ':00:00'        ####
            print(file_path)
            for i in range(start_lon, end_lon + 1, 1):
                lon_temp = Left_lon_temp + step * i
                for j in range(start_lat,end_lat + 1, 1):
                    lat_temp = Down_lat_temp + step * j
                    temp_val = data.iloc[j,i]
                    if math.isnan(temp_val):
                        continue
                    if float(lat_temp) > 41 and float(lat_temp) < 43:
                        North += temp_val
                    elif float(lon_temp) < 115:
                        West += temp_val
                    elif float(lat_temp) < 38.85:
                        South += temp_val
                    elif float(lon_temp) < 117:
                        local += temp_val
                    else:
                        East += temp_val

            res = pd.concat([pd.Series(date), pd.Series(local),pd.Series(North),pd.Series(West),pd.Series(South),pd.Series(East)], axis=1, ignore_index=False)
            result.append(res)
    header = ['Time', 'Local', 'North', 'West', 'South', 'East']
    result = pd.concat(result, ignore_index=False)
    result.to_csv('/Volumes/HK/NUIST/CWT/CWT/results_2013.csv', sep=',', header = header, index=False)

def value_to_radians():
    pwd = '/Volumes/HK/NUIST/CWT/CWT/time_ECMWF.txt'
    f = open('/Volumes/HK/NUIST/CWT/CWT/Res_DOY_ECMWF.csv', 'a')
    data = open(pwd, 'r')
    data.readline()
    while 1:
        line = data.readline()
        if not line:
            break
        year = int(line.strip().split(' ')[0].split('/')[0])
        month_original = int(line.strip().split(' ')[0].split('/')[1])
        line_split = line.strip()
        time1 = datetime.strptime(line_split, "%Y/%m/%d %H:%M:%S")
        DOY_original = int(time1.timetuple().tm_yday)
        DOW_original = int(time1.strftime("%w")) + 1
        if (year % 4 == 0 and year % 100 != 0) or year % 400 == 0:
            DOY = (math.pi / 182.5 * (DOY_original - 183.5))
        else:
            DOY = (math.pi / 183 * (DOY_original - 184))
        DOW = math.pi / 3.5 * (DOW_original - 4.5)
        Month = math.pi / 6 * (month_original-7)
        f.write('%10s,%10s,%10s,%10s,%10s,%10s,%10s\n' % (year, DOY_original, DOW_original,month_original, DOY, DOW, Month))
    f.close()
    data.close()




     #print (math.pi / 6 * (value-7))        ####Month
    # print(math.pi / 11.5 * (value - 12.5))          ####Hour
    #print(math.pi / 3.5 * (value - 4.5))            ####DOW
     #print(math.pi / 182.5 * (value - 183.5))        ####DOY 非闰年
    #print(math.pi / 183 * (value - 184))  ####DOY 闰年

#-----------------------------------------------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    value_to_radians()


