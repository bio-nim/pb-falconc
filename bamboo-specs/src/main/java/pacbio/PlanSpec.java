package pacbio;

import com.atlassian.bamboo.specs.api.BambooSpec;
import com.atlassian.bamboo.specs.api.builders.AtlassianModule;
import com.atlassian.bamboo.specs.api.builders.BambooKey;
import com.atlassian.bamboo.specs.api.builders.BambooOid;
import com.atlassian.bamboo.specs.api.builders.permission.PermissionType;
import com.atlassian.bamboo.specs.api.builders.permission.Permissions;
import com.atlassian.bamboo.specs.api.builders.permission.PlanPermissions;
import com.atlassian.bamboo.specs.api.builders.plan.Job;
import com.atlassian.bamboo.specs.api.builders.plan.Plan;
import com.atlassian.bamboo.specs.api.builders.plan.PlanIdentifier;
import com.atlassian.bamboo.specs.api.builders.plan.Stage;
import com.atlassian.bamboo.specs.api.builders.plan.branches.BranchCleanup;
import com.atlassian.bamboo.specs.api.builders.plan.branches.PlanBranchManagement;
import com.atlassian.bamboo.specs.api.builders.plan.configuration.ConcurrentBuilds;
import com.atlassian.bamboo.specs.api.builders.project.Project;
import com.atlassian.bamboo.specs.api.builders.requirement.Requirement;
import com.atlassian.bamboo.specs.api.builders.task.AnyTask;
import com.atlassian.bamboo.specs.builders.task.CheckoutItem;
import com.atlassian.bamboo.specs.builders.task.ScriptTask;
import com.atlassian.bamboo.specs.builders.task.TestParserTask;
import com.atlassian.bamboo.specs.builders.task.VcsCheckoutTask;
import com.atlassian.bamboo.specs.builders.trigger.BitbucketServerTrigger;
import com.atlassian.bamboo.specs.model.task.ScriptTaskProperties;
import com.atlassian.bamboo.specs.model.task.TestParserTaskProperties;
import com.atlassian.bamboo.specs.util.BambooServer;
import com.atlassian.bamboo.specs.util.MapBuilder;

@BambooSpec
public class PlanSpec {
    
    public Plan plan() {
        final Plan plan = new Plan(new Project()
             
                 .key(new BambooKey("SAT"))
                .name("SMRT Analysis Tools (SAT)"),
            "falconc",
            new BambooKey("FALCONC"))
          
            .description("Plan created from Bamboo Java Specs,  modify http://bitbucket.pacificbiosciences.com:7990/projects/SAT/repos/falconc/browse project to update the plan.")
            .pluginConfigurations(new ConcurrentBuilds()
                    .useSystemWideDefault(false)
                    .maximumNumberOfConcurrentBuilds(4))
            .stages(new Stage("Default Stage")
                    .jobs(new Job("Default Job",
                            new BambooKey("JOB1"))
                            .tasks(new VcsCheckoutTask()
                                    .description("Checkout Default Repository")
                                    .checkoutItems(new CheckoutItem().defaultRepository()),
                                new ScriptTask()
                                    .description("build")
                                    .location(ScriptTaskProperties.Location.FILE)
                                    .fileFromPath("bamboo_build.sh"))
                            .finalTasks(new TestParserTask(TestParserTaskProperties.TestType.JUNIT)
                                    .description("junitparsertry")
                                    .resultDirectories("tests/falconc_unit_tests_junit.xml"),
                                new AnyTask(new AtlassianModule("com.xpandit.plugins.xray-bamboo:cucumberFeaturesImportTask"))
                                    .description("import falconc unit tests into X-ray")
                                    .enabled(false)
                                    .configuration(new MapBuilder()
                                            .put("format", "")
                                            .put("importServer", "2")
                                            .put("inputInfoSwitcher", "")
                                            .put("importFilePath", "tests/xray_summary.json")
                                            .build()))
                            .requirements(new Requirement("system.os")
                                    .matchValue("linux")
                                    .matchType(Requirement.MatchType.EQUALS))))
            .linkedRepositories("falconc")
            
            .triggers(new BitbucketServerTrigger())
            .planBranchManagement(new PlanBranchManagement()
                    .createForPullRequest()
                    .delete(new BranchCleanup()
                        .whenRemovedFromRepositoryAfterDays(7)
                        .whenInactiveInRepositoryAfterDays(30))
                    .notificationForCommitters())
            .forceStopHungBuilds();
        return plan;
    }
    
    public PlanPermissions planPermission() {
        final PlanPermissions planPermission = new PlanPermissions(new PlanIdentifier("SAT", "FALCONC"))
            .permissions(new Permissions()
                    .userPermissions("cdunn", PermissionType.BUILD, PermissionType.CLONE, PermissionType.ADMIN, PermissionType.VIEW, PermissionType.EDIT)
                    .userPermissions("zkronenberg", PermissionType.BUILD, PermissionType.CLONE, PermissionType.ADMIN, PermissionType.VIEW, PermissionType.EDIT)
                    .userPermissions("bli", PermissionType.EDIT, PermissionType.BUILD, PermissionType.CLONE, PermissionType.VIEW, PermissionType.ADMIN));
        return planPermission;
    }
    
    public static void main(String... argv) {
        //By default credentials are read from the '.credentials' file.
        BambooServer bambooServer = new BambooServer("http://bamboo.pacificbiosciences.com:8085");
        final PlanSpec planSpec = new PlanSpec();
        
        final Plan plan = planSpec.plan();
        bambooServer.publish(plan);
        
        final PlanPermissions planPermission = planSpec.planPermission();
        bambooServer.publish(planPermission);
    }
}
