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
    /// flag.xaml の相互作用ロジック
    /// </summary>
    public partial class toggleButton : UserControl
    {
        public toggleButton()
        {
            InitializeComponent();
        }
        public string Text
        {
            set
            {
                this.button1.Content = value;
            }
            get
            {
                return (string)this.button1.Content;
            }
        }
        public bool value
        {
            get
            {
                return (bool)this.button1.IsChecked;
            }
        }
        public static implicit operator bool(Mothra.UI.toggleButton b)
        {
            return b.value;
        }
    }
}
