import gdb

# 計數器，用來區分不同階段的檔案
count = 0

def my_handler(event):
    global count
    # 檢查是否停在我們想要的檔案和行號
    frame = gdb.selected_frame()
    sal = frame.find_sal()
    
    if sal.symtab and "CUFLU_FluidSolver_MHM.cu" in sal.symtab.filename and sal.line == 508:
        filename = f"flux_dump_{count}.txt"
        print(f"\n[AutoDump] 偵測到斷點，正在寫入 {filename}...")
        
        gdb.execute(f"set logging file {filename}")
        gdb.execute("set logging overwrite on")
        gdb.execute("set logging on")
        
        # 執行取值
        gdb.execute("x/10648gf (@global real *)&g_Flux_Half_1PG[2][12][0]")
        
        gdb.execute("set logging enabled off")
        print(f"[AutoDump] {filename} 儲存成功。")
        count += 1
        
        # 存完後自動繼續執行
        gdb.execute("continue")

# 註冊事件監聽器（當程式停止時觸發）
gdb.events.stop.connect(my_handler)

# 設定初始環境
gdb.execute("set pagination off")
gdb.execute("break CUFLU_FluidSolver_MHM")
gdb.execute("run")
gdb.execute("break CUFLU_FluidSolver_MHM.cu:508")
gdb.execute("continue")