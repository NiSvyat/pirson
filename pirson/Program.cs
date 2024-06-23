using System;

class pirson{
  static void Main() 
  {
    int N=100000; //количество отсчетов
    double sigma_2 = 0.25;  //объявление константы
    double[] x1 = new double[N]; //массив для последовательности чисел х1
    double[] x2 = new double[N]; //массив для последовательности чисел х2
    double[] ksi1 = new double[N]; //массив для последовательности чисел ksi1
    double[] ksi2 = new double[N]; //массив для последовательности чисел ksi2
    double[] t = new double[N]; //массив для коррелированной последовательности тау
    double U1=0, U2=0, ln=0, pi=0; //объявление и обнуление вспомогательных переменных
    double sum, sum_x1, sum_x2, sum_x1_2, sum_x2_2, sum_tx1, sum_tx2, sum_t_2, sum_x1x2; //объявление переменных для подсчета сумм
    double r_2; //объявление переменной r в квадрате
      
    double ksi1p, ksi2p; //объявление переменных "ksi1 предыдущее" и "ksi2 предыдущее"
    double runner1, runner2, result; //объявление переменных для расчета коэффициента корреляции по Пирсону
    
    Random r = new Random(); //создание переменной для генерации рандомных чисел
    
    //цикл для моделирования последовательности с коэффициентов корреляции от 0.5 до 0.95 с шагом 0.05 
    for(double p1 = 0.5; p1 < 1; p1+= 0.05)
    {
        r_2 = p1; //присвоение переменной r в квадрате значения коэффициента корреляции
        
        //обнуление сумм
        sum = 0;
        sum_x1 = 0;
        sum_x2 = 0;
        sum_x1_2 = 0;
        sum_x2_2 = 0;
        sum_tx1 = 0;
        sum_tx2 = 0;
        sum_t_2 = 0;
        sum_x1x2 = 0;
         
        //моделирование значений x1[0] и x2[0]
        U1 = r.NextDouble();
        U2 = r.NextDouble();
        pi = 2*Math.PI*U2;
        ln = Math.Sqrt(-2*Math.Log(U1));
        x1[0] = ln*Math.Cos(pi);
        x2[0] = ln*Math.Sin(pi);
        
        //моделирование первых значений последовательностей ksi1 и ksi2
        ksi1[0] = sigma_2*Math.Pow((Math.Sqrt(1-r_2)*x1[0]), 2.0);
        ksi2[0] = sigma_2*Math.Pow((Math.Sqrt(1-r_2)*x2[0]), 2.0);

        //присвоение переменным ksi1 предыдущее и ksi2 предыдущее текущих значений ksi1 и ksi2
        ksi1p = ksi1[0];
        ksi2p = ksi2[0];

        //моделирование первого отсчета (первого занчения коррелированной последовательности)
        t[0] = ksi1[0] + ksi2[0];
    
        //добавление к суммам первых значений последовательностей
        sum_tx1 += t[0] * ksi1[0];
        sum_tx2 += t[0] * ksi2[0];
        sum_x1x2 += ksi1[0] * ksi2[0];
        sum +=t[0]; 
  
        sum_x1 += ksi1[0];
        sum_x2 += ksi2[0];
              
        sum_t_2+=t[0] * t[0];
        sum_x1_2+=(ksi1[0]*ksi1[0]);
        sum_x2_2+=(ksi2[0]*ksi2[0]);
        
        //цикл для моделирования следующих значений последовательностей
        for(int i=1; i<N; i++)
        {
            U1 = r.NextDouble();
            U2 = r.NextDouble();
            pi = 2*Math.PI*U2;
            ln = Math.Sqrt(-2*Math.Log(U1));
            x1[i] = ln*Math.Cos(pi);
            x2[i] = ln*Math.Sin(pi);
              
            ksi1[i] = sigma_2 * Math.Pow( Math.Sqrt( 1.0 - r_2) * x1[i] + ( Math.Sqrt(r_2) * Math.Sqrt(ksi1p)), 2.0);
            ksi1[i] = sigma_2 * Math.Pow( Math.Sqrt( 1.0 - r_2) * x2[i] + ( Math.Sqrt(r_2) * Math.Sqrt(ksi2p)), 2.0);
              
            ksi1p = ksi1[i];
            ksi2p = ksi2[i];
              
            t[i]=ksi1[i]+ksi2[i];
              
            sum_tx1 += t[i] * ksi1[i];
            sum_tx2 += t[i] * ksi2[i];
            sum_x1x2 += ksi1[i] * ksi2[i];
            sum +=t[i];
              
            sum_x1 += ksi1[i];
            sum_x2 += ksi2[i];
              
            sum_t_2+=t[i] * t[i];
            sum_x1_2+=(ksi1[i]*ksi1[i]);
            sum_x2_2+=(ksi2[i]*ksi2[i]);
        }

        //расчет и вывод зависимости тау от кси 1
        runner1 = N * sum_tx1 - sum * sum_x1;
        runner2 = Math.Sqrt( (N * sum_t_2 - Math.Pow(sum, 2) ) * (N * sum_x1_2 - Math.Pow(sum_x1, 2) )   );
        result = runner1/runner2;
        Console.WriteLine("Зависимость Тау от Кси 1 = " + result + " при коэффициенте корреляции p = " + p1);
          
        //расчет и вывод зависимости тау от кси 2
        runner1 = N * sum_tx2 - sum * sum_x2;
        runner2 = Math.Sqrt( (N * sum_t_2 - Math.Pow(sum, 2) ) * (N * sum_x2_2 - Math.Pow(sum_x2, 2) )   );
        result = runner1/runner2;
        Console.WriteLine("Зависимость Тау от Кси 2 = " + result + " при коэффициенте корреляции p = " + p1);
        
        //расчет и вывод зависимости кси 1 от кси 2
        runner1 = ( N * sum_x1x2 ) - (sum_x1 * sum_x2);
        runner2 = Math.Sqrt( (N * sum_x1_2 - Math.Pow(sum_x1, 2) ) * (N * sum_x2_2 - Math.Pow(sum_x2, 2) )   );
        result = runner1/runner2;
        Console.WriteLine("Зависимость Кси 1 от Кси 2 = " + result + " при коэффициенте корреляции p = " + p1);
        }
    }
}