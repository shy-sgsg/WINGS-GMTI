/*
 * PipeRW.h
 *
 *  Created on: Feb 23, 2021
 *      Author: raco
 */

#ifndef PIPERW_H_
#define PIPERW_H_

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <stdint.h>
#include <iostream>

#pragma pack(1)

typedef enum eGmtiWorkMode
{
    Mode_SAR_01 = 1,
    Mode_SAR_03 = 2,
    Mode_SAR_05 = 3,
    Mode_SAR_1  = 4,
    Mode_SAR_3  = 5,
    Mode_GMTI   = 7,     // GMTI mode value

    Mode_SINGLEPULSE = 8, //单脉冲扫描

    Mode_INIT   = 100,    //初始上电
}eGmtiWorkMode;

#pragma pack()

class PipeRW {
public:
	PipeRW();
	virtual ~PipeRW();

    bool CreatePipe(const char* fifoFileName, int flag);
	void ClosePipe();

    ssize_t  WriteData(const char* pData, int nLen);
	int  ReadData(char* pBuf, int nBufSize);

private:
	char m_fifoFileName[256];
	int m_fd;
	int m_pid;
};

#endif /* PIPERW_H_ */
