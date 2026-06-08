#ifndef PIPESTRUDEF_H
#define PIPESTRUDEF_H

#include <sys/types.h>
#include <stdint.h>
#include <stdio.h>

#pragma pack(push, 1)
// 发送参数及工作命令结构体
typedef struct ModeSwitchCmd
{
    uint32_t head;
    uint8_t algoType;
    uint32_t dataLen;
    uint32_t workMode;
    double CenterFreq;
    double SamplingRate;
    double PulseWidth;
    double BandWidth;
    double prf;
    uint32_t Kr_Sign;
    double SampleDelay;
    uint32_t updatePointNum;
    uint32_t reserved0;
    uint32_t LookSide;
    double LookDownAngle;      //下视角
    double SquintAngle;        //斜视角
    double Theta_bw;           //方位波束宽度
    double alt_scene;          //成像区域高度（地面参考海拔）
    float Rmin;                //最近斜距
    uint64_t rfFreq;           //工作载频，和中心频率重复
    double Kr;                 //调频率
    double ADsamplingLen;      //AD采样波门宽度
    uint32_t channelSpace;     //通道间隔
    float velocity;            //速度
    uint16_t swath;            //幅宽
    double antennaAz;          //天线方位角
    double antennaEl;          //天线俯仰角
    double targetLon;          //目标区域经度
    double targetLat;          //目标区域纬度
    double targetAlt;          //目标区域高度
    uint16_t targetRange;      //目标斜距
    uint8_t reserved1;
    uint32_t tail;
}ModeSwitchCmd;


typedef struct ResultHeader
{
    uint16_t head;
    uint32_t msgLen;
    uint16_t msgAddr;
    uint16_t msgType;
    uint16_t msgCount;
    uint16_t srcId;
    uint16_t dstId;
    uint8_t cmdType;
    uint8_t cmdCount;
    uint16_t height;
    uint16_t width;
    uint16_t availFlag;

    int16_t roll;
    int16_t heading;
    int16_t pitch;
    int32_t navLon;
    int32_t navLat;
    int16_t navHeight;
    int16_t velNorth;
    int16_t velUp;
    int16_t velEast;
    uint8_t hour; //
    uint8_t minute; //
    uint8_t second; //
    uint8_t millisec; //
    uint8_t posReserve[48];

    int16_t hLeftTop;
    int16_t hLeftDown;
    int16_t hRightDown;
    int16_t hRightTop;
    int16_t hCenter;
    int32_t lonLeftTop;
    int32_t lonLeftDown;
    int32_t lonRightDown;
    int32_t lonRightTop;
    int32_t lonCenter;
    int32_t latLeftTop;
    int32_t latLeftDown;
    int32_t latRightDown;
    int32_t latRightTop;
    int32_t latCenter;
    uint16_t rLeftTop;
    uint16_t rLeftDown;
    uint16_t rRightDown;
    uint16_t rRightTop;
    uint16_t rCenter;
    uint8_t reserve1;
    uint8_t pixelPitch; //
    uint16_t LookDownAngle;
    uint16_t SquintAngle;
    uint8_t LookSide;
    uint8_t reserve2[4];
    uint8_t checksum;
    uint16_t targetNum;
}ResultHeader;

typedef struct EchoFrameHead
{
    uint8_t res1[9];
    uint32_t PRTlen;
    uint8_t res2[243];
}EchoFrameHead;

#pragma pack(pop)

#endif // PIPESTRUDEF_H
