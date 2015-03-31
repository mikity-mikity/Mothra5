using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace Mothra.UI
{

    /// <summary>
    /// variable.xaml の相互作用ロジック
    /// </summary>
    public partial class rangeSetter : UserControl
    {
        private string key;
        Action<rangeSetter,int, double, double> func;
        private void setKey(string _key) {
            key = _key;
            this.Label1.Content = key + ":lb";
            this.Label2.Content = key + ":ub";
        }
        public rangeSetter()
        {
            InitializeComponent();
        }
        public rangeSetter(string key)
        {
            InitializeComponent();
            setKey(key);
        }
        public void setFunc(Action<rangeSetter,int, double, double> _func)
        {
            func = _func;
        }
        public void setMeasured(double min,double max){
            this.measured1.Content = min.ToString("g4");
            this.measured2.Content = max.ToString("g4");
        }
        private void CheckBox_Checked(object sender, RoutedEventArgs e)
        {

            if (ubEnabel.IsChecked == true)
            {
                double lb, ub;
                bool isLb=true, isUb=true;
                try
                {
                    lb = double.Parse(this.lbValue.Text);
                    
                }
                catch
                {
                    lb = 0;
                    isLb = false;
                }
                try
                {
                    ub = double.Parse(this.ubValue.Text);
                }
                catch
                {
                    ub = 0;
                    isUb = false;
                }

                if (func != null)
                {
                    if (isUb && isLb)
                        func(this,2, lb, ub);
                    if (isUb && (!isLb))
                        func(this,2, 0d, ub);
                    if (!(isUb) && isLb)
                        func(this,0, lb, 0d);
                    if (!(isUb) && !(isLb))
                        func(this,0, 0d, 0d);
                }
            }
            else
            {
                double lb;
                bool isLb = true, isUb = false;
                try
                {
                    lb = double.Parse(this.lbValue.Text);

                }
                catch
                {
                    lb = 0;
                    isLb = false;
                }

                if (func != null)
                {
                    if (!(isUb) && isLb)
                        func(this,0, lb, 0d);
                    if (!(isUb) && !(isLb))
                        func(this,0, 0d, 0d);
                }

            }
        }

        private void CheckBox_Unchecked(object sender, RoutedEventArgs e)
        {
            if (ubEnabel.IsChecked == true)
            {
                double ub;
                bool isLb = false, isUb = true;
                try
                {
                    ub = double.Parse(this.ubValue.Text);

                }
                catch
                {
                    ub = 0;
                    isUb = false;
                }

                if (func != null)
                {
                    if (isUb && (!isLb))
                        func(this,2, 0d, ub);
                    if (!(isUb) && !(isLb))
                        func(this,0, 0d, 0d);
                }
            } else
            {
                if (func != null)
                {
                    func(this,0, 0d, 0d);
                }
            }

        }

        private void CheckBox_Unchecked_1(object sender, RoutedEventArgs e)
        {
            if (lbEnabel.IsChecked == true)
            {
                double lb, ub;
                bool isLb = true, isUb = true;
                try
                {
                    lb = double.Parse(this.lbValue.Text);

                }
                catch
                {
                    lb = 0;
                    isLb = false;
                }
                try
                {
                    ub = double.Parse(this.ubValue.Text);
                }
                catch
                {
                    ub = 0;
                    isUb = false;
                }

                if (func != null)
                {
                    if (isUb && isLb)
                        func(this,2, lb, ub);
                    if (isUb && (!isLb))
                        func(this,2, 0d, ub);
                    if (!(isUb) && isLb)
                        func(this,0, lb, 0d);
                    if (!(isUb) && !(isLb))
                        func(this,0, 0d, 0d);
                }
            }
            else
            {
                double ub;
                bool isUb = true, isLb = false;
                try
                {
                    ub = double.Parse(this.ubValue.Text);

                }
                catch
                {
                    ub = 0;
                    isUb = false;
                }

                if (func != null)
                {
                    if (!(isLb) && isUb)
                        func(this,2, 0d, ub);
                    if (!(isUb) && !(isLb))
                        func(this,0, 0d, 0d);
                }

            }

        }

        private void CheckBox_Checked_1(object sender, RoutedEventArgs e)
        {
            if (lbEnabel.IsChecked == true)
            {
                double lb;
                bool isLb = true, isUb = false;
                try
                {
                    lb = double.Parse(this.lbValue.Text);
                }
                catch
                {
                    lb = 0;
                    isLb = false;
                }

                if (func != null)
                {
                    if ((!isUb) && (isLb))
                        func(this,0, lb, 0d);
                    if (!(isUb) && !(isLb))
                        func(this,0, 0d, 0d);
                }
            }

            else
            {
                if (func != null)
                {
                    func(this,0, 0d, 0d);
                }
            }
        }
        public void update()
        {
            double lb, ub;
            bool isLb = true, isUb = true;
            if (lbEnabel.IsChecked == true)
            {
                try
                {
                    lb = double.Parse(this.lbValue.Text);

                }
                catch
                {
                    lb = 0;
                    isLb = false;
                }
            }
            else
            {
                isLb = false;
                lb = 0;
            }
            if (ubEnabel.IsChecked == true)
            {
                try
                {
                    ub = double.Parse(this.ubValue.Text);

                }
                catch
                {
                    ub = 0;
                    isUb = false;
                }
            }
            else
            {
                isUb = false;
                ub = 0;
            }

            if (func != null)
            {
                if (isUb && isLb)
                    func(this,2, lb, ub);
                if (isUb && (!isLb))
                    func(this,2, 0d, ub);
                if (!(isUb) && isLb)
                    func(this,0, lb, 0d);
                if (!(isUb) && !(isLb))
                    func(this,0, 0d, 0d);
            }

        }
    }
}

