using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SPSA
{
    public class BCD
    {
        public static double FuncDescend(int loca, int scope, double[] x)
        {
            return loca switch
            {
                0 => Math.Sqrt((Math.Pow(x[loca + 1], 2) + Math.Pow(x[loca + 2], 2)) / 3),
                1 when scope - loca - 1 >= 0 => Math.Sqrt((Math.Pow(x[loca - 1], 2) + 2 * Math.Pow(x[loca + 1], 2) + Math.Pow(x[loca + 2], 2)) / 6),
                >= 2 when scope - loca - 1 >= 0 => Math.Sqrt((Math.Pow(x[loca - 2], 2) + 2 * Math.Pow(x[loca - 1], 2) + 2 * Math.Pow(x[loca + 1], 2) + Math.Pow(x[loca + 2], 2)) / 9),
                1 when scope - loca - 1 == -1 => Math.Sqrt((Math.Pow(x[loca + 1], 2) + Math.Pow(x[loca - 1], 2)) / 3),
                >= 2 when scope - loca - 1 == -1 => Math.Sqrt((Math.Pow(x[loca - 2], 2) + 2 * Math.Pow(x[loca - 1], 2) + Math.Pow(x[loca + 1], 2)) / 6),
                >= 2 when scope - loca - 1 == -2 => Math.Sqrt((Math.Pow(x[loca - 2], 2) + Math.Pow(x[loca - 1], 2)) / 3),
                _ => -1
            };
        }
        public double FuncXi(int loca, int scope, double[] x) 
        {
            double fi = 0;
            for (int i = 0; i < scope; i++)
            {
                fi += -2 * x[i] * x[i + 1] - 2 * x[i] * x[i + 2] - 2 * x[i + 2] * x[i + 1];
            }
            return fi;
        }
    }
}
