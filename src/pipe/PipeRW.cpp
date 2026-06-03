/*
 * PipeRW.cpp
 *
 *  Created on: Feb 23, 2021
 *      Author: raco
 */

#include "PipeRW.h"
#include "PipeStruDef.h"

PipeRW::PipeRW() {
	// TODO Auto-generated constructor stub
	m_pid = -1;
	m_fd = -1;
}

PipeRW::~PipeRW() {
	// TODO Auto-generated destructor stub
}

bool PipeRW::CreatePipe(const char* fifoFileName, int flag)
{
    memset(m_fifoFileName, 0, sizeof(m_fifoFileName));
	memcpy(m_fifoFileName, fifoFileName, strlen(fifoFileName));

	m_pid = mkfifo(m_fifoFileName, 0777);
    if(m_pid == -1)
	{
        if(errno == EEXIST)
        {
            std::cout << "fifo channel " << m_fifoFileName << " exists" << std::endl;
        }
        else
        {
        	perror("create fifo channel failed.");
        	return false;
        }
	}

    m_fd = open(m_fifoFileName, flag);
    if(m_fd == -1)
	{
		perror("open fifo failed.");
		return false;
	}
    if (fcntl(m_fd, F_SETPIPE_SZ, 1 * 1024 * 1024) == -1) {
        perror("fcntl");
        // 错误处理
    }

	return true;
}

void PipeRW::ClosePipe()
{
    if(m_fd == -1)
	{
		return;
	}
	close(m_fd);
}

ssize_t  PipeRW::WriteData(const char* pData, int nLen)
{
    if(m_fd == -1)
	{
		return 0;
	}

	return write(m_fd, pData, nLen);
}

int PipeRW::ReadData(char* pBuf, int nBufSize)
{
    if(m_fd == -1)
	{
		return 0;
	}

	return read(m_fd, pBuf, nBufSize);
}

