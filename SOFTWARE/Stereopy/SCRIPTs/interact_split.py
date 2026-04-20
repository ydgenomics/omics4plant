# 启动交互式空间分布图
# 在弹出的界面中，侧边栏工具栏里有一个“Lasso Select”工具
st.plt.interact_spatial(data, res_key='leiden', width=600, height=600)
st.plt.interact_spatial(data, res_key='None', width=600, height=600)