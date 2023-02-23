# Technical Indicators

技术指标基础算法,只依赖于C++标准库,类似麦语言语义.

### 用法示例:

#### MACD

    //DIF:EMA(CLOSE,SHORT)-EMA(CLOSE,LONG);
    //DEA:EMA(DIF,MID);
    //MACD:(DIF-DEA)*2,COLORSTICK;

    std::vector<double> CLOSE = ...;
    unsigned int SHORT = ...;
    unsigned int LONG = ...;
    unsigned int MID = ...;

    std::vector<double> DIF = vu_Sub(vu_EMA(CLOSE,SHORT),vu_EMA(CLOSE,LONG));
    std::vector<double> DEA = vu_EMA(DIF,MID);
    std::vector<double> macd = vu_Mul(2,vu_Sub(DIF,DEA));


