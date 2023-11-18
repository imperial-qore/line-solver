package jline.lib.thirdparty;

import jline.util.Matrix;

import java.util.HashMap;
import java.util.Map;

public class SMCSolver {
    public static Map<String, Matrix> QBD_CR(Matrix A0, Matrix A1, Matrix A2, Integer MaxNumIt_, Integer Verbose_, String Mode_, Integer RAPComp_){
        String Mode = "Shift";
        int MaxNumIt = 50;
        boolean Verbose = false;
        boolean RAPComp = false;

        if(MaxNumIt_!=null){
            MaxNumIt = MaxNumIt_;
        }

        if(Mode_!=null){
            if (Mode_.equals("Shift")||Mode_.equals("Basic")) {
                Mode = Mode_;
            }else {
                throw new RuntimeException("QBD_LR mode not recognized");
            }
        }

        if(Verbose_!=null){
            if(Verbose_==1){
                Verbose = true;
            }
        }

        if(RAPComp_!=null){
            if(RAPComp_==1){
                RAPComp = true;
            }
        }
        Matrix A1_diag = new Matrix(0,0,0);
        Matrix.extractDiag(A1,A1_diag);
        boolean continues = true;
        double lamb = Matrix.negative(A1_diag).elementMax();
        int m = A1.numRows;
        if(!RAPComp){
            continues = false;
            if(A1_diag.elementSum()<0){
                continues =true;
                A0.scale(1/lamb);
                A1.scale(1/lamb);
                A1 = A1.add(1,Matrix.eye(m));
                A2.scale(1/lamb);
            }
            QBD_ParsePara(A0,A1,A2);
        }else{
            QBD_ParsePara(A0,A1,A2);
            A0.scale(1/lamb);
            A1.scale(1/lamb);
            A1 = A1.add(1,Matrix.eye(m));
            A2.scale(1/lamb);
        }

        Map<String,Matrix> result = QBD_EG(A0,A1,A2,Verbose);
        Matrix G = result.get("G");
        if(G.length()==0){
            return result;
        }
        Matrix theta = stat(A0.add(1,A1).add(1,A2));
        double drift = theta.mult(A0.sumRows()).get(0)-theta.mult(A2.sumRows()).get(0);
        Matrix A2old = A2.clone();
        Matrix uT = Matrix.scale_mult(Matrix.ones(1,m),m);
        Matrix A0old = A0.clone();
        if(Mode.equals("Shift")){
            if(drift<0){
                A2 = A2.add(-1,Matrix.ones(m,1).mult(theta.mult(A2)));
                A1 = A1.add(1,Matrix.ones(m,1).mult(theta.mult(A0)));
            }else {
                A0 = A0.add(-1,A0.sumRows().mult(uT));
                A1 = A1.add(1,A2.sumRows().mult(uT));
            }
        }
        Matrix A = A1.clone();
        Matrix B = A2.clone();
        Matrix C = A0.clone();

        Matrix Ahat = A.clone();
        double check =1;
        int numit = 0;

        while (check>Math.pow(10,-14)&&numit<MaxNumIt){
            Matrix Atemp = Matrix.eye(m).add(-1,A).inv();
            Matrix BAtemp = B.mult(Atemp);
            Atemp = C.mult(Atemp);
            Ahat = A.add(1,BAtemp.mult(C));
            A = A.add(1,BAtemp.mult(A)).add(1,Atemp.mult(B));
            B = BAtemp.mult(B);
            C = Atemp.mult(C);
            numit = numit+1;
            check = Math.min(Matrix.inf_norm(B),Matrix.inf_norm(C));
            if(Verbose){
                System.out.println("Check after "+numit+"iterations: "+ check );
            }
        }
        if(numit == MaxNumIt && check>Math.pow(10,-14)){
            System.out.println("Maximum Number of Iterations reached");
        }

        G = Matrix.eye(m).add(-1,Ahat).inv().mult(A0);

        if(Mode.equals ("Shift")){
            if(drift<0){
                A1 = A1.add(-1,Matrix.ones(m,1).mult(theta).mult(A0));
                A2 = A2old.clone();
            }else {
                G = G.add(1,Matrix.ones(m,1).mult(uT));
                A1 = A1.add(-1,A2.sumRows().mult(uT));
                A0 = A0old.clone();
            }
        }
        if (Verbose){
            double res_norm = Matrix.inf_norm(G.add(-1,A0.add(-1,A1.add(1,A2.mult(G)).mult(G))));
            System.out.println("Final Residual Error for G: "+ res_norm);
        }

        Matrix R = A2.mult(Matrix.eye(m).add(-1,A1.add(1,A2.mult(G))).inv());
        if(Verbose){
            double res_norm = Matrix.inf_norm(R.add(-1,A2).add(-1,R.mult(A1.add(1,R.mult(A0)))));
            System.out.println("Final Residual Error for R:" + res_norm);
        }

        Matrix U = A1.add(1,R.mult(A0));
        if(Verbose){
            double res_norm = Matrix.inf_norm(U.add(-1,A1).add(-1,A2.mult(Matrix.eye(m).add(-1,U).inv()).mult(A0)));
            System.out.println("Final Residual Error for U " + res_norm);
        }
        if(continues){
            U = Matrix.scale_mult(U.add(-1,Matrix.eye(m)),lamb);
        }

        result = new HashMap<>();
        result.put("G",G);
        result.put("R",R);
        result.put("U",U);
        return result;
    }

    public static Map<String,Matrix> QBD_LR(Matrix A0, Matrix A1, Matrix A2, Integer MaxNumIt_, Integer Verbose_, String Mode_, Integer RAPComp_){
        String Mode = "Shift";
        int MaxNumIt = 50;
        boolean Verbose = false;
        boolean RAPComp = false;

        if(MaxNumIt_!=null){
            MaxNumIt = MaxNumIt_;
        }

        if(Mode_!=null){
            if (Mode_.equals("Shift")||Mode_.equals("Basic")) {
                Mode = Mode_;
            }else {
                throw new RuntimeException("QBD_LR mode not recognized");
            }
        }

        if(Verbose_!=null){
            if(Verbose_==1){
                Verbose = true;
            }
        }

        if(RAPComp_!=null){
            if(RAPComp_==1){
                RAPComp = true;
            }
        }
        Matrix A1_diag = new Matrix(0,0,0);
        Matrix.extractDiag(A1,A1_diag);
        boolean continues = true;
        double lamb = Matrix.negative(A1_diag).elementMax();
        int m = A1.numRows;
        if(!RAPComp){
            continues = false;
            if(A1_diag.elementSum()<0){
                continues =true;
                A0.scale(1/lamb);
                A1.scale(1/lamb);
                A1 = A1.add(1,Matrix.eye(m));
                A2.scale(1/lamb);
            }
            QBD_ParsePara(A0,A1,A2);
        }else{
            QBD_ParsePara(A0,A1,A2);
            A0.scale(1/lamb);
            A1.scale(1/lamb);
            A1 = A1.add(1,Matrix.eye(m));
            A2.scale(1/lamb);
        }
        Map<String,Matrix> result = QBD_EG(A0,A1,A2,Verbose);
        Matrix G = result.get("G");
        if(G.length()==0){
            return result;
        }
        Matrix theta = stat(A0.add(1,A1).add(1,A2));
        double drift = theta.mult(A0.sumRows()).get(0)-theta.mult(A2.sumRows()).get(0);
        Matrix A2old = A2.clone();
        Matrix uT = Matrix.scale_mult(Matrix.ones(1,m),m);
        Matrix A0old = A0.clone();
        if(Mode.equals("Shift")){
            if(drift<0){
                A2 = A2.add(-1,Matrix.ones(m,1).mult(theta.mult(A2)));
                A1 = A1.add(1,Matrix.ones(m,1).mult(theta.mult(A0)));
            }else {
                A0 = A0.add(-1,A0.sumRows().mult(uT));
                A1 = A1.add(1,A2.sumRows().mult(uT));
            }
        }

        Matrix B2 = Matrix.eye(m).add(-1,A1).inv();
        Matrix B0 = B2.mult(A2);
        B2 =  B2.mult(A0);
        G = B2.clone();
        Matrix PI = B0.clone();
        double check =1 ;
        int numit =0;
        while (check>Math.pow(10,-14)&& numit<MaxNumIt){
            Matrix A1star = B2.mult(B0).add(1,B0.mult(B2));
            Matrix A0star = B0.mult(B0);
            Matrix A2star = B2.mult(B2);
            B0 = Matrix.eye(m).add(-1,A1star).inv();
            B2 = B0.mult(A2star);
            B0 = B0.mult(A0star);
            G = G.add(1, PI.mult(B2));
            PI = PI.mult(B0);
            check = Math.min(Matrix.inf_norm(B0),Matrix.inf_norm(B2));
            numit = numit+1;
            if (Verbose){
                System.out.println("Check after "+numit+" iterations: "+check);
            }
        }

        if(numit == MaxNumIt && check>Math.pow(10,-14)){
            System.out.println("Maximum Number of Iterations reached");
        }

        if(Mode.equals ("Shift")){
            if(drift<0){
                A1 = A1.add(-1,Matrix.ones(m,1).mult(theta).mult(A0));
                A2 = A2old.clone();
            }else {
                G = G.add(1,Matrix.ones(m,1).mult(uT));
                A1 = A1.add(-1,A2.sumRows().mult(uT));
                A0 = A0old.clone();
            }
        }
        if (Verbose){
            double res_norm = Matrix.inf_norm(G.add(-1,A0.add(-1,A1.add(1,A2.mult(G)).mult(G))));
            System.out.println("Final Residual Error for G: "+ res_norm);
        }

        Matrix R = A2.mult(Matrix.eye(m).add(-1,A1.add(1,A2.mult(G))).inv());
        if(Verbose){
            double res_norm = Matrix.inf_norm(R.add(-1,A2).add(-1,R.mult(A1.add(1,R.mult(A0)))));
            System.out.println("Final Residual Error for R:" + res_norm);
        }

        Matrix U = A1.add(1,R.mult(A0));
        if(Verbose){
            double res_norm = Matrix.inf_norm(U.add(-1,A1).add(-1,A2.mult(Matrix.eye(m).add(-1,U).inv()).mult(A0)));
            System.out.println("Final Residual Error for U " + res_norm);
        }
        if(continues){
            U = Matrix.scale_mult(U.add(-1,Matrix.eye(m)),lamb);
        }

        result = new HashMap<>();
        result.put("G",G);
        result.put("R",R);
        result.put("U",U);
        return result;
    }

    public static Map<String,Matrix> QBD_NI(Matrix A0, Matrix A1, Matrix A2, Integer MaxNumIt_, Integer Verbose_, String Mode_, Integer RAPComp_){
        String Mode = "Sylvest";
        int MaxNumIt = 50;
        boolean Verbose = false;
        boolean RAPComp = false;

        if(MaxNumIt_!=null){
            MaxNumIt = MaxNumIt_;
        }

        if(Mode_!=null){
            if (Mode_.equals("Sylvest")||Mode_.equals("Estimat")||Mode_.equals("DirectSum")) {
                Mode = Mode_;
            }else {
                throw new RuntimeException("QBD_LR mode not recognized");
            }
        }

        if(Verbose_!=null){
            if(Verbose_==1){
                Verbose = true;
            }
        }

        if(RAPComp_!=null){
            if(RAPComp_==1){
                RAPComp = true;
            }
        }
        Matrix A1_diag = new Matrix(0,0,0);
        Matrix.extractDiag(A1,A1_diag);
        boolean continues = true;
        double lamb = Matrix.negative(A1_diag).elementMax();
        int m = A1.numRows;
        if(!RAPComp){
            continues = false;
            if(A1_diag.elementSum()<0){
                continues =true;
                A0.scale(1/lamb);
                A1.scale(1/lamb);
                A1 = A1.add(1,Matrix.eye(m));
                A2.scale(1/lamb);
            }
            QBD_ParsePara(A0,A1,A2);
        }else{
            QBD_ParsePara(A0,A1,A2);
            A0.scale(1/lamb);
            A1.scale(1/lamb);
            A1 = A1.add(1,Matrix.eye(m));
            A2.scale(1/lamb);
        }
        Map<String,Matrix> result = QBD_EG(A0,A1,A2,Verbose);
        Matrix G = result.get("G");
        if(G.length()==0){
            return result;
        }
        Matrix R = new Matrix(m,m,m*m);
        double check =1;
        int numit =0;
        while (check>Math.pow(10,-12)&&numit<MaxNumIt){
            numit = numit+1;
            Matrix YK;
            if(numit==1){
                YK = A2.mult(Matrix.eye(m).add(-1,A1).inv());
            }else {
                if(Mode.equals("Estimat")){
                    Matrix FRK = A2.add(1,R.mult(A1.add(-1,Matrix.eye(m)).add(1,R.mult(A0))));
                    Matrix ZK = FRK.mult(Matrix.eye(m).add(-1,A1)).inv();
                    YK = FRK.add(1,FRK.add(1,ZK.mult(A1)).add(1,R.mult(ZK).add(1,ZK.mult(R)).mult(A0)));
                }else {
                    Matrix D = Matrix.scale_mult(A2.add(1,R.mult(A1.add(-1,Matrix.eye(m)).add(1,R.mult(A0)))),-1);
                    Matrix C = A1.add(1,R.mult(A0)).add(-1,Matrix.eye(m));
                    if(Mode.equals("Sylvest")){
                        YK = QBD_NI_Sylvest(A0.transpose(),R.transpose(),C.transpose(),D.transpose());
                    }else {
                        Matrix D_reshape = D.clone();
                        D_reshape.reshape(m*m,1);
                        YK = A0.transpose().krons(R).add(1,C.transpose().krons(Matrix.eye(m)).inv().mult(D_reshape));
                        YK.reshape(m,m);
                    }
                }
            }
            R = R.add(1,YK);
            check = Matrix.inf_norm(YK);
            if (Verbose){
                System.out.println("Check after "+numit+" iterations: "+check);
            }
        }
        if(numit == MaxNumIt && check>Math.pow(10,-12)){
            System.out.println("Maximum Number of Iterations reached");
        }

        if (Verbose){
            double res_norm = Matrix.inf_norm(G.add(-1,A0.add(-1,A1.add(1,A2.mult(G)).mult(G))));
            System.out.println("Final Residual Error for G: "+ res_norm);
        }

        R = A2.mult(Matrix.eye(m).add(-1,A1.add(1,A2.mult(G))).inv());
        if(Verbose){
            double res_norm = Matrix.inf_norm(R.add(-1,A2).add(-1,R.mult(A1.add(1,R.mult(A0)))));
            System.out.println("Final Residual Error for R:" + res_norm);
        }

        Matrix U = A1.add(1,R.mult(A0));
        if(Verbose){
            double res_norm = Matrix.inf_norm(U.add(-1,A1).add(-1,A2.mult(Matrix.eye(m).add(-1,U).inv()).mult(A0)));
            System.out.println("Final Residual Error for U " + res_norm);
        }
        if(continues){
            U = Matrix.scale_mult(U.add(-1,Matrix.eye(m)),lamb);
        }

        result = new HashMap<>();
        result.put("G",G);
        result.put("R",R);
        result.put("U",U);
        return result;
    }

    public static Map<String,Matrix> QBD_FI(Matrix A0, Matrix A1, Matrix A2, Integer MaxNumIt_, Integer Verbose_, String Mode_, Matrix StartValue_, Integer RAPComp_){
        String Mode = "U-Based";
        int MaxNumIt = 10000;
        boolean Verbose = false;
        boolean RAPComp = false;
        int m = A1.numRows;
        Matrix StartValue = new Matrix(m,m,m*m);

        if(StartValue_!=null){
            StartValue = StartValue_.clone();
        }

        if(MaxNumIt_!=null){
            MaxNumIt = MaxNumIt_;
        }

        if(Mode_!=null){
            if (Mode_.equals("Traditional ")||Mode_.equals("Natural")||Mode_.equals("U-Based")||Mode_.equals("ShiftTraditional")||Mode_.equals("ShiftNatural")||Mode_.equals("ShiftU-Based")) {
                Mode = Mode_;
            }else {
                throw new RuntimeException("QBD_LR mode not recognized");
            }
        }

        if(Verbose_!=null){
            if(Verbose_==1){
                Verbose = true;
            }
        }

        if(RAPComp_!=null){
            if(RAPComp_==1){
                RAPComp = true;
            }
        }
        Matrix A1_diag = new Matrix(0,0,0);
        Matrix.extractDiag(A1,A1_diag);
        boolean continues = true;
        double lamb = Matrix.negative(A1_diag).elementMax();

        if(!RAPComp){
            continues = false;
            if(A1_diag.elementSum()<0){
                continues =true;
                A0.scale(1/lamb);
                A1.scale(1/lamb);
                A1 = A1.add(1,Matrix.eye(m));
                A2.scale(1/lamb);
            }
            QBD_ParsePara(A0,A1,A2);
        }else{
            QBD_ParsePara(A0,A1,A2);
            A0.scale(1/lamb);
            A1.scale(1/lamb);
            A1 = A1.add(1,Matrix.eye(m));
            A2.scale(1/lamb);
        }
        Map<String,Matrix> result = QBD_EG(A0,A1,A2,Verbose);
        Matrix G = result.get("G");
        if(G.length()==0){
            return result;
        }
        int numit =0;
        double check =1;
        G = StartValue.clone();
        Matrix theta = stat(A0.add(1,A1).add(1,A2));
        double drift = theta.mult(A0.sumRows()).get(0)-theta.mult(A2.sumRows()).get(0);
        Matrix A2old = A2.clone();
        Matrix uT = Matrix.scale_mult(Matrix.ones(1,m),m);
        Matrix A0old = A0.clone();
        if(Mode.contains("Shift")){
            if(drift<0){
                A2 = A2.add(-1,Matrix.ones(m,1).mult(theta.mult(A2)));
                A1 = A1.add(1,Matrix.ones(m,1).mult(theta.mult(A0)));
            }else {
                A0 = A0.add(-1,A0.sumRows().mult(uT));
                A1 = A1.add(1,A2.sumRows().mult(uT));
            }
        }

        if(Mode.contains("Natural")){
            while (check>Math.pow(10,-14)&&numit<MaxNumIt){
                Matrix Gold = G.clone();
                G = A2.mult(G).add(1,A1).mult(G).add(1,A0);
                check = Matrix.inf_norm(G.add(-1,Gold));
                numit = numit+1;
                if (Verbose){
                    System.out.println("Check after "+numit+" iterations: "+check);
                }
            }
        }

        if(Mode.contains("Traditional")){
            Matrix invA1 = Matrix.eye(m).add(-1,A1).inv();
            while (check>Math.pow(10,-14)&&numit<MaxNumIt){
                Matrix Gold = G.clone();
                G = invA1.mult(A0.add(1,A2.mult(Matrix.pow(G,2))));
                check = Matrix.inf_norm(G.add(-1,Gold));
                numit = numit+1;
                if (Verbose){
                    System.out.println("Check after "+numit+" iterations: "+check);
                }
            }
        }

        if(Mode.contains("U-Based")){
            while (check>Math.pow(10,-14)&&numit<MaxNumIt){
                Matrix Gold = G.clone();
                G = Matrix.eye(m).add(-1,A1).add(-1,A2.mult(G)).inv().mult(A0);
                check = Matrix.inf_norm(G.add(-1,Gold));
                numit = numit+1;
                if (Verbose){
                    System.out.println("Check after "+numit+" iterations: "+check);
                }
            }
        }

        if(numit == MaxNumIt && check>Math.pow(10,-12)){
            System.out.println("Maximum Number of Iterations reached");
        }

        if(Mode.contains ("Shift")){
            if(drift<0){
                A1 = A1.add(-1,Matrix.ones(m,1).mult(theta).mult(A0));
                A2 = A2old.clone();
            }else {
                G = G.add(1,Matrix.ones(m,1).mult(uT));
                A1 = A1.add(-1,A2.sumRows().mult(uT));
                A0 = A0old.clone();
            }
        }
        if (Verbose){
            double res_norm = Matrix.inf_norm(G.add(-1,A0.add(-1,A1.add(1,A2.mult(G)).mult(G))));
            System.out.println("Final Residual Error for G: "+ res_norm);
        }

        Matrix R = A2.mult(Matrix.eye(m).add(-1,A1.add(1,A2.mult(G))).inv());
        if(Verbose){
            double res_norm = Matrix.inf_norm(R.add(-1,A2).add(-1,R.mult(A1.add(1,R.mult(A0)))));
            System.out.println("Final Residual Error for R:" + res_norm);
        }

        Matrix U = A1.add(1,R.mult(A0));
        if(Verbose){
            double res_norm = Matrix.inf_norm(U.add(-1,A1).add(-1,A2.mult(Matrix.eye(m).add(-1,U).inv()).mult(A0)));
            System.out.println("Final Residual Error for U " + res_norm);
        }
        if(continues){
            U = Matrix.scale_mult(U.add(-1,Matrix.eye(m)),lamb);
        }

        result = new HashMap<>();
        result.put("G",G);
        result.put("R",R);
        result.put("U",U);
        return result;
    }

    public static Map<String,Matrix> QBD_IS(Matrix A0, Matrix A1, Matrix A2, Integer MaxNumIt_, Integer Verbose_, String Mode_, Integer RAPComp_){
        String Mode = "Schur";
        int MaxNumIt = 50;
        boolean Verbose = false;
        boolean RAPComp = false;
        int m = A1.numRows;

        if(MaxNumIt_!=null){
            MaxNumIt = MaxNumIt_;
        }

        if(Mode_!=null){
            if (Mode_.equals("MSignStandard")||Mode_.equals("MSignBalzer")||Mode_.equals("Schur")) {
                Mode = Mode_;
            }else {
                throw new RuntimeException("QBD_LR mode not recognized");
            }
        }

        if(Verbose_!=null){
            if(Verbose_==1){
                Verbose = true;
            }
        }

        if(RAPComp_!=null){
            if(RAPComp_==1){
                RAPComp = true;
            }
        }
        Matrix A1_diag = new Matrix(0,0,0);
        Matrix.extractDiag(A1,A1_diag);
        boolean continues = true;
        double lamb = Matrix.negative(A1_diag).elementMax();

        if(!RAPComp){
            continues = false;
            if(A1_diag.elementSum()<0){
                continues =true;
                A0.scale(1/lamb);
                A1.scale(1/lamb);
                A1 = A1.add(1,Matrix.eye(m));
                A2.scale(1/lamb);
            }
            QBD_ParsePara(A0,A1,A2);
        }else{
            QBD_ParsePara(A0,A1,A2);
            A0.scale(1/lamb);
            A1.scale(1/lamb);
            A1 = A1.add(1,Matrix.eye(m));
            A2.scale(1/lamb);
        }
        Map<String,Matrix> result = QBD_EG(A0,A1,A2,Verbose);
        Matrix G = result.get("G");
        if(G.length()==0){
            return result;
        }

        double epsilon = Math.pow(10,-12);
        int f = 2;

        Matrix theta = stat(A0.add(1,A1).add(1,A2));
        double drift = theta.mult(A0.sumRows()).get(0)-theta.mult(A2.sumRows()).get(0);
        Map<Integer,Matrix> F = new HashMap<>();
        F.put(1,Matrix.scale_mult(A0,-1));
        F.put(2,Matrix.eye(m).add(-1,A1));
        F.put(3,Matrix.scale_mult(A2,-1));
        Map<Integer,Matrix> H = new HashMap<>();
        for(int i=0;i<=f;i++){
            H.put(i+1,new Matrix(m,m));
        }
        for(int i=0;i<=f;i++){
            Matrix temp1 = new Matrix(1,2,2);
            temp1.set(0,0,1);
            temp1.set(0,1,-1);
            Matrix temp2 = new Matrix(1,2,2);
            temp2.set(0,0,1);
            temp2.set(0,1,1);
        }

        //todo convolution,orth,schur,hesse

        if (Verbose){
            double res_norm = Matrix.inf_norm(G.add(-1,A0.add(-1,A1.add(1,A2.mult(G)).mult(G))));
            System.out.println("Final Residual Error for G: "+ res_norm);
        }

        Matrix R = A2.mult(Matrix.eye(m).add(-1,A1.add(1,A2.mult(G))).inv());
        if(Verbose){
            double res_norm = Matrix.inf_norm(R.add(-1,A2).add(-1,R.mult(A1.add(1,R.mult(A0)))));
            System.out.println("Final Residual Error for R:" + res_norm);
        }

        Matrix U = A1.add(1,R.mult(A0));
        if(Verbose){
            double res_norm = Matrix.inf_norm(U.add(-1,A1).add(-1,A2.mult(Matrix.eye(m).add(-1,U).inv()).mult(A0)));
            System.out.println("Final Residual Error for U " + res_norm);
        }
        if(continues){
            U = Matrix.scale_mult(U.add(-1,Matrix.eye(m)),lamb);
        }

        result = new HashMap<>();
        result.put("G",G);
        result.put("R",R);
        result.put("U",U);
        return result;
    }

    public static void QBD_ParsePara(Matrix A0, Matrix A1, Matrix A2){
        if(A0.numCols!=A0.numRows){
            throw new RuntimeException("MATLAB:QBD_ParsePara:InvalidInput, A0 is not a square matrix");
        }
        if(A1.numCols!=A1.numRows){
            throw new RuntimeException("MATLAB:QBD_ParsePara:InvalidInput, A1 is not a square matrix");
        }
        if(A2.numCols!=A2.numRows){
            throw new RuntimeException("MATLAB:QBD_ParsePara:InvalidInput, A2 is not a square matrix");
        }
        if(A0.numCols!=A1.numCols){
            throw new RuntimeException("MATLAB:QBD_ParsePara:InvalidInput, The matrices A0 and A1 do not have the same dimension");
        }
        if(A0.numCols!=A2.numCols){
            throw new RuntimeException("MATLAB:QBD_ParsePara:InvalidInput, The matrices A0 and A2 do not have the same dimension");
        }
        if(A0.elementMin()<=-Math.pow(-10,-14)){
            throw new RuntimeException("MATLAB:QBD_ParsePara:InvalidInput, The matrix A0 contains negative data");
        }
        if(A1.elementMin()<=-Math.pow(-10,-14)){
            throw new RuntimeException("MATLAB:QBD_ParsePara:InvalidInput, The matrix A1 contains negative data");
        }
        if(A2.elementMin()<=-Math.pow(-10,-14)){
            throw new RuntimeException("MATLAB:QBD_ParsePara:InvalidInput, The matrix A2 contains negative data");
        }
        if(A0.add(1,A1).add(1,A2).sumRows().elementMax()>1+Math.pow(10,-14)){
            throw new RuntimeException("MATLAB:QBD_ParsePara:InvalidInput, The matrix A0+A1+A2 has to be (sub)stochastic");
        }
    }

    public static Map<String,Matrix> QBD_EG(Matrix A0, Matrix A1, Matrix A2,Boolean optVerbose){
        //todo define matrix pow
        Matrix G = new Matrix(0,0,0);
        Matrix R = new Matrix(0,0,0);
        Matrix U = new Matrix(0,0,0);
        int m = A1.numRows;
        Matrix theta = stat(A0.add(1,A1).add(1,A2));
        double drift= theta.safe_mult(A0.sumRows()).get(0) - theta.safe_mult(A2.sumRows()).get(0);
        if(drift >0){
            if(A0.rank()==1){
                // matlab rank
                int non_zero_row = 0;
                Matrix row_sum = A0.sumRows();
                for(int i=0;i<A0.numRows;i++){
                    if(row_sum.get(i)>0){
                        non_zero_row = i;
                        break;
                    }
                }
                Matrix beta = Matrix.scale_mult(Matrix.extractRows(A0,non_zero_row,non_zero_row+1,null),1/A0.sumRows(non_zero_row));
                G = Matrix.ones(m,1).mult(beta);
                R = A2.mult(Matrix.eye(m).add(-1,A1.add(1,A2.mult(G))).inv());
            }else if(A2.rank()==1){
                double eta = QBD_CAUDAL(A0, A1, A2);
                R = A2.mult(Matrix.eye(m).add(-1,A1).add(-eta,A0).inv());
                G = Matrix.eye(m).add(-1,A1.add(1,R.mult(A0))).inv().mult(A0);
            }
        }else if(drift<0){
            if(A2.rank()==1){
                Matrix alpha = A2.mult(Matrix.ones(m,1));
                R = Matrix.scale_mult(alpha.mult(theta),theta.mult(alpha).get(0));
                G = Matrix.eye(m).add(-1,A1.add(1,R.mult(A0))).inv().mult(A0);
            }else if(A0.rank()==0){
                Matrix A0hat = Matrix.diag(theta.element_power(-1).toArray1D()).mult(A2.transpose()).mult(Matrix.diag(theta.toArray1D()));
                Matrix A1hat = Matrix.diag(theta.element_power(-1).toArray1D()).mult(A1.transpose()).mult(Matrix.diag(theta.toArray1D()));
                Matrix A2hat = Matrix.diag(theta.element_power(-1).toArray1D()).mult(A0.transpose()).mult(Matrix.diag(theta.toArray1D()));
                double etahat = QBD_CAUDAL(A0hat,A1hat,A2hat);
                G = Matrix.diag(theta.element_power(-1).toArray1D()).mult(A2hat.mult(Matrix.eye(m).add(-1,A1hat).add(-etahat,A0hat).inv()).transpose()).mult(Matrix.diag(theta.toArray1D()));
                R = A2.mult(Matrix.eye(m).add(-1,A1.add(1,A2.mult(G))).inv());
            }
        }

        if(R.length()>0){
            U = A1.add(1,R.mult(A0));
        }

        if(optVerbose){
            if(G.length()>0){
                double resnorm = Matrix.inf_norm(G.add(-1,A0).add(-1,A1.add(1,A2.mult(G)).mult(G)));
                System.out.println("Final Residual Error for G: "+resnorm);
            }
            if(R.length()>0){
                double resnorm = Matrix.inf_norm(R.add(-1,A2).add(-1,R.mult(A1.add(1,R.mult(A0)))));
                System.out.println("Final Residual Error for R: "+resnorm);
            }
            if(U.length()>0){
                double resnorm = Matrix.inf_norm(U.add(-1,A1).add(-1,A2).mult(Matrix.eye(m).add(-1,U).inv()).mult(A0));
                System.out.println("Final Residual Error for U: "+resnorm);
            }
        }
        Map<String,Matrix> result = new HashMap<>();
        result.put("G",G);
        result.put("R",R);
        result.put("U",U);
        return result;
    }

    public static double QBD_CAUDAL (Matrix A0, Matrix A1, Matrix A2, Boolean Dual){
        if(Dual){
            Matrix A2old = A2;
            A2 = A0;
            A0 = A2old;
        }

        double eta_min = 0;
        double eta_max = 1;
        double eta = 0.5;
        while(eta_max-eta_min>Math.pow(10,15)){
            double new_eta = A2.add(1,Matrix.scale_mult(A1,eta)).add(1,Matrix.scale_mult(A0,Math.pow(eta,2))).eigenvalue().elementMax();
            if(new_eta>eta){
                eta_min = eta;
            }else {
                eta_max = eta;
            }
            eta = (eta_max+eta_min)/2;
        }
        return eta;
    }

    public static double QBD_CAUDAL(Matrix A0, Matrix A1, Matrix A2){
        return QBD_CAUDAL(A0,A1,A2,false);
    }

    public static Matrix stat(Matrix A){
        int S = A.numRows;
        Matrix e = Matrix.ones(S,1);
        Matrix B = Matrix.concatColumns(A.add(-1,Matrix.eye(S)),e,null);
        Matrix y = Matrix.concatColumns(new Matrix(1,S),new Matrix(1),null);
        return  B.right_matrix_divide(y);
    }




    public static Matrix QBD_NI_Sylvest(Matrix A, Matrix B, Matrix C, Matrix D){
        throw new RuntimeException("QBD_NI_Sylvest has not been implemented");
    }
}
