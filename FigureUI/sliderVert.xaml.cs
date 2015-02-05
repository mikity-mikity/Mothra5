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
    public partial class sliderVert : UserControl
    {
        public System.Windows.Controls.Slider getSlider
        {
            get { return this.slider1; }
        }
        public System.Windows.Controls.Label getLabel
        {
            get { return this.label1; }
        }
        public string text;
        private Func<double, double> _convert;
        public int originalValue
        {
            get
            {
                return (int)this.slider1.Value;
            }
        }
        public double value
        {
            get
            {
                return _convert(this.slider1.Value);
            }
        }
        public sliderVert()
        {
            InitializeComponent();
        }
        public sliderVert(int min, int step, int max, int val, string _text)
        {
            InitializeComponent();
            this.slider1.Minimum = min;
            this.slider1.Maximum = max;
            this.slider1.Value = val;
            for (int i = min; i <= max; i += step)
            {
                this.slider1.Ticks.Add(i);
            }
            _convert = (v) => { return v; };
            this.text = _text;
            this.update();
        }
        public void update()
        {
            if (Math.Abs(this.value) < 0.01)
            {
                this.label1.Content = this.text + ":  " + this.value.ToString("0.000e00");
            }
            else
            {
                this.label1.Content = this.text + ":  " + this.value.ToString("0.000");
            }
        }
        public Action<bool> fixChanged;
        public Func<double, double> Converter
        {
            set
            {
                _convert = value;
                update();
            }
        }


    }
}

