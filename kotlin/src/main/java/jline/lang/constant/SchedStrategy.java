package jline.lang.constant;

/**
 *  Constants for specifying scheduling strategies at stations
 */
public enum SchedStrategy {
    INF,
    FCFS,
    LCFS,
    LCFSPR,
    SIRO,
    SJF,
    LJF,
    PS,
    DPS,
    GPS,
    SEPT,
    LEPT,
    HOL,
    FORK,
    EXT,
    REF;

    public static SchedStrategy fromLINEString(String string) {
        switch(string) {
            case "inf":
                return INF;
            case "fcfs":
                return FCFS;
            case "lcfs":
                return LCFS;
            case "lcfspr":
                 return LCFSPR;
            case "siro":
                return SIRO;
            case "sjf":
                return SJF;
            case "ljf":
                return LJF;
            case "ps":
                return PS;
            case "dps":
                return DPS;
            case "gps":
                return GPS;
            case "sept":
                return SEPT;
            case "lept":
                return LEPT;
            case "hol":
                return HOL;
            case "fork":
                return FORK;
            case "ext":
                return EXT;
            case "ref":
                return REF;
            default:
                throw new RuntimeException("Unable to return a SchedStrategy, check string and try again.");
        }
    }


    public static String toText(SchedStrategy scheduling) {
        switch (scheduling){
            case INF:
                return "inf";
            case FCFS:
                return "fcfs";
            case LCFS:
                return "lcfs";
            case LCFSPR:
                return "lcfspr";
            case SIRO:
                return "siro";
            case PS:
                return "ps";
            case DPS:
                return "dps";
            case GPS:
                return "gps";
            case SEPT:
                return "sept";
            case LEPT:
                return "lept";
            case HOL:
                return "hol";
            case FORK:
                return "fork";
            case EXT:
                return "ext";
            case REF:
                return "ref";
            default:
                return "";
        }
    }

    public static String toFeature(SchedStrategy scheduling) {
        switch (scheduling){
            case INF:
                return "SchedStrategy_INF";
            case FCFS:
                return "SchedStrategy_FCFS";
            case LCFS:
                return "SchedStrategy_LCFS";
            case LCFSPR:
                return "SchedStrategy_LCFSPR";
            case SIRO:
                return "SchedStrategy_SIRO";
            case SJF:
                return "SchedStrategy_SJF";
            case LJF:
                return "SchedStrategy_LJF";
            case PS:
                return "SchedStrategy_PS";
            case DPS:
                return "SchedStrategy_DPS";
            case GPS:
                return "SchedStrategy_GPS";
            case SEPT:
                return "SchedStrategy_SEPT";
            case LEPT:
                return "SchedStrategy_LEPT";
            case HOL:
                return "SchedStrategy_HOL";
            case EXT:
                return "SchedStrategy_EXT";
            default:
                return "";
        }
    }

    public static int toID(SchedStrategy scheduling) {
        switch (scheduling){
            case INF:
                return 0;
            case FCFS:
                return 1;
            case LCFS:
                return 2;
            case SIRO:
                return 3;
            case SJF:
                return 4;
            case LJF:
                return 5;
            case PS:
                return 6;
            case DPS:
                return 7;
            case GPS:
                return 8;
            case SEPT:
                return 9;
            case LEPT:
                return 10;
            case HOL:
                return 11;
            case FORK:
                return 12;
            case EXT:
                return 13;
            case REF:
                return 14;
            case LCFSPR:
                return 15;
            default:
                return -1;
        }
    }
}