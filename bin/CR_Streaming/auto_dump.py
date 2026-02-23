import gdb
import os

FIELDS = {
    "CR_E": 12,
    "CR_F1": 11,
    "CR_F2": 10,
    "CR_F3": 9,
    "ADV_SIGMA": 8,
}

# ================= 配置區域：在這邊增加你要抓取的點 =================
# 格式: (檔案名, 行號, 抓取表達式, 元素數量, 存檔別名)
MHM_CU = "CUFLU_FluidSolver_MHM.cu"
DUMP_TARGETS = [
    (MHM_CU, 508, "&g_Flux_Half_1PG[2]", 22**3, "Flux_Half_2"),
    (MHM_CU, 518, "&g_Flu_Array_In[P]",  22**3, "Flu_Array_In"),
    (MHM_CU, 544, "&g_PriVar_Half_1PG",  22**3, "PriVar_Half"), # Stride 20, 0 offset
    # (MHM_CU, 593, "&g_FC_Var_1PG",       18**3, "FC_Var_1PG"),
    (MHM_CU, 621, "&g_FC_Flux_1PG[2]",   22**3, "FC_Flux_1PG_2"), # Stride 18
    # 你可以隨時增加更多點...
]

# DUMP_TARGETS = [
#     (MHM_CU, 447, "&g_Flu_Array_In[P]", 22**3, "01_Flu_Array_In_start"),
#     (MHM_CU, 481, "&g_Flu_Array_In[P]", 22**3, "02_Flu_Array_In_UpdateOpacity_F"),
#     (MHM_CU, 508, "&g_Flu_Array_In[P]", 22**3, "03_Flu_Array_In_UpdateStreaming_H"),
#     (MHM_CU, 529, "&g_Flu_Array_In[P]", 22**3, "04_Flu_Array_In_Hydro_RiemannPredict_H"),
#     (MHM_CU, 655, "&g_Flu_Array_In[P]", 22**3, "05_Flu_Array_In_Hydro_FullStepUpdate_F"),
# ]


DUMP_DIR = "dump/raw"
# =================================================================

class SmartDumpManager:
    def __init__(self, targets):
        self.targets = targets
        self.hit_counts = {i: 0 for i in range(len(targets))}
        
        # 建立資料夾
        if not os.path.exists(DUMP_DIR):
            os.makedirs(DUMP_DIR)
            
        # 註冊斷點
        for i, (filename, line, expr, size, tag) in enumerate(self.targets):
            gdb.execute(f"break {filename}:{line}")
            print(f"[Init] 已對 {filename}:{line} 設定自動監控")

    def handle_stop(self, event):
        # 取得目前停下的位置
        frame = gdb.selected_frame()
        sal = frame.find_sal()
        
        if not sal.symtab:
            return

        current_file = os.path.basename(sal.symtab.filename)
        current_line = sal.line

        # 檢查是否符合我們的目標清單
        for i, (filename, line, expr, size, tag) in enumerate(self.targets):
            if current_file == filename and current_line == line:
                self.execute_dump(i)
            elif current_file == filename and current_line == 380: # MHM 第一行指令
                gdb.execute("continue")

    def execute_dump(self, index):
        filename, line, expr, size, tag = self.targets[index]
        step = self.hit_counts[index] // 8
        patch = self.hit_counts[index] % 8
        
        print(f"[SmartDump] 命中: {tag} (Step {step}, Patch {patch}) -> 正在背景取值...", end='', flush=True)
        
        # --- 改良部分：使用 to_string=True 隱藏輸出 ---
        try:
            for field, shift in FIELDS.items():
                cmd = f"x/{size}gf (@global real *){expr}[{shift}][0]"
                raw_output = gdb.execute(cmd, to_string=True)
                
                # 直接用 Python 把字串存入檔案
                output_path = os.path.join(DUMP_DIR, f"{tag}_{field}_step{step}_patch{patch}.txt")
                with open(output_path, "w") as f:
                    f.write(raw_output)
            
            print(f" [完成]")
        except Exception as e:
            print(f" [失敗: {e}]")
        # --------------------------------------------
        
        self.hit_counts[index] += 1
        
        # 自動繼續執行
        gdb.execute("continue")

# 初始化
if __name__ == "__main__":
    manager = SmartDumpManager(DUMP_TARGETS)
    gdb.events.stop.connect(manager.handle_stop)
    print("[SmartDump] 系統啟動。")
