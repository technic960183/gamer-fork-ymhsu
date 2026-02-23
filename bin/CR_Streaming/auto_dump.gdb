# 1. 基本設定
set pagination off

# 2. 設定第一個斷點 (進入 Function)
break CUFLU_FluidSolver_MHM
run

# 3. 設定第二個斷點 (目標取值行號)
# 我們為這個斷點編號 2，並綁定自動指令
break CUFLU_FluidSolver_MHM.cu:508
commands 2
  echo --- 正在執行自動取值作業 ---\n
  set logging file flux_dump.txt
  set logging overwrite on
  set logging on
  
  # 執行你成功的取值指令
  x/10648gf (@global real *)&g_Flux_Half_1PG[2][12][0]
  
  set logging off
  echo --- 取值完畢，結果已存入 flux_dump.txt ---\n
  
  # 如果你想要取完值後讓程式自動往下跑，可以取消下面 continue 的註釋
  # continue
end

# 4. 讓程式繼續跑到 508 行
continue